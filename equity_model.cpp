#include<iostream>
#include"equity_model.h"
#include<ql/quantlib.hpp>

using namespace std;
using namespace QuantLib;
 
void heston_model_european_options() { 
	// option paramters
	Real strike_price = 130.0;
	auto payoff = ext::make_shared<PlainVanillaPayoff>(Option::Call, strike_price);
	
	//option data
	Date maturity_date(15, January, 2016);
	Real spot_price = 127.62;
	Real volatility = 0.20;
	Real dividend_rate = 0.0163;
	Option::Type option_type = Option::Call;

	Real risk_free_rate = 0.001;
	DayCounter day_count = Actual365Fixed();
	Calendar calendar = UnitedStates(UnitedStates::NERC);

	Date calculation_date(8, May, 2015);
	Settings::instance().evaluationDate() = calculation_date;

	auto exercise = ext::make_shared<EuropeanExercise>(maturity_date);
	VanillaOption european_option(payoff, exercise);


	Real v0 = volatility * volatility;
	Real kappa = 0.1, theta = v0, sigma = 0.1, rho = -0.75;
	Handle<Quote> spot_handle(ext::make_shared<SimpleQuote>(spot_price));
	Handle<YieldTermStructure> flat_ts(ext::make_shared<FlatForward>(calculation_date, risk_free_rate, day_count));
	Handle<YieldTermStructure> dividend_yield(ext::make_shared<FlatForward>(calculation_date, dividend_rate, day_count));
	auto heston_process = ext::make_shared<HestonProcess>(flat_ts, dividend_yield, spot_handle, v0, kappa, theta, sigma, rho);

	auto engine = ext::make_shared<AnalyticHestonEngine>(ext::make_shared<HestonModel>(heston_process), 0.01, 1000);
	european_option.setPricingEngine(engine);
	Real h_price = european_option.NPV();
	cout << "The Heston model price is " << h_price << endl;


	Handle<BlackVolTermStructure> flat_vol_ts(ext::make_shared<BlackConstantVol>(calculation_date, calendar, volatility, day_count));
	auto bsm_process = ext::make_shared<BlackScholesMertonProcess>(spot_handle, dividend_yield, flat_ts, flat_vol_ts);
	european_option.setPricingEngine(ext::make_shared<AnalyticEuropeanEngine>(bsm_process));
	Real bs_price = european_option.NPV();
	cout << "The Black-Scholes-Merton model price is " << bs_price << endl;

}
// binomial tree pricing
Real binominal_price(VanillaOption & option, const ext::shared_ptr<BlackScholesMertonProcess> &bsm_process, const int &steps) {
	auto engine = ext::make_shared<BinomialVanillaEngine<CoxRossRubinstein>> (bsm_process,steps);
	option.setPricingEngine(engine); 
	return option.NPV();
	
}
void valuing_european_american_options() {
	Date maturity_date(15, January, 2016);
	Real spot_price = 127.62, strike_price = 130, volatility = 0.20, dividend_rate = 0.0163;
	Option::Type option_type = Option::Call;

	Real risk_free_rate = 0.001;
	DayCounter day_count = Actual365Fixed();
	Calendar calendar = UnitedStates(UnitedStates::Settlement);
	Date calculation_date(8, May, 2015);

	Settings::instance().evaluationDate() = calculation_date;
	auto payoff = ext::make_shared<PlainVanillaPayoff>(option_type, strike_price);
	auto exercise = ext::make_shared<EuropeanExercise>(maturity_date);
	VanillaOption european_option(payoff, exercise);

	Handle<Quote> spot_handle(ext::make_shared<SimpleQuote>(spot_price));

	Handle<YieldTermStructure> flat_ts(ext::make_shared<FlatForward>(calculation_date, risk_free_rate, day_count));
	Handle<YieldTermStructure> dividend_yield(ext::make_shared<FlatForward>(calculation_date, dividend_rate, day_count));
	Handle<BlackVolTermStructure> flat_vol_ts(ext::make_shared<BlackConstantVol>(calculation_date, calendar, volatility, day_count));
	auto bsm_process = ext::make_shared<BlackScholesMertonProcess>(spot_handle, dividend_yield, flat_ts, flat_vol_ts);
	european_option.setPricingEngine(ext::make_shared<AnalyticEuropeanEngine>(bsm_process));
	Real bs_price = european_option.NPV();
	cout << "the theoretical price is : " << bs_price << endl;

	vector<Real> prices;
	for (int i = 2; i < 200; i++) {
		prices.push_back(binominal_price(european_option, bsm_process, i));
	}
//	for (auto& price : prices) { cout << price << endl; }
	
	//American Options

	
	Date settlement = calculation_date;
	auto am_exercise = ext::make_shared<AmericanExercise>(settlement, maturity_date);
	VanillaOption american_option(payoff, am_exercise);
	prices.clear();
	for (int i = 2; i < 200; i++) {
		prices.push_back(binominal_price(american_option, bsm_process, i));
	}
	//for (auto& price : prices) { cout << price << endl; }

}

void black_formula_pricing() {
	Calendar calendar = UnitedStates(UnitedStates::SOFR);
	BusinessDayConvention business_convention = ModifiedFollowing;
	Natural settlement_days = 0;
	DayCounter day_count = ActualActual(ActualActual::Bond);

	//option on Treasury Futures Contract
	Real interest_rate = 0.00105;
	Date calc_date(1, December, 2015);
	auto yield_curve = ext::make_shared<FlatForward>(calc_date, interest_rate, day_count, Compounded);

	Settings::instance().evaluationDate() = calc_date;
	Date option_maturity_date(24, December, 2015);
	Real strike = 119;
	Real spot = 126.953;
	Real volatility = 11.567 / 100;
	Option::Type flavor = Option::Call;
	Real discount = yield_curve->discount(option_maturity_date);
	auto strikepayoff = ext::make_shared<PlainVanillaPayoff>(flavor, strike);
	Time t = yield_curve->dayCounter().yearFraction(calc_date, option_maturity_date);
	Real stddev = volatility * sqrt(t);

	BlackCalculator black(strikepayoff, spot, stddev, discount); 
	cout << "Option Price: " << black.value() << endl;
	cout << "Delta: " << black.delta(spot) << endl;
	cout << "Gamma: " << black.gamma(spot) << endl;
	cout << "Theta: " << black.theta(spot, t) << endl;
	cout << "Vega: " << black.vega(t) << endl;
	cout << "Rho: " << black.rho(t) << endl; 

	//Natural Gas Futures Option
	cout << endl << "Natural Gas Futures Option " << endl;

	interest_rate = 0.0015;
	calc_date = Date(23, September, 2015);
	yield_curve = ext::make_shared<FlatForward>(calc_date, interest_rate, day_count, Compounded); 
	Settings::instance().evaluationDate() = calc_date;
	Real T = 96.12 / 365;
	strike = 3.5;
	spot = 2.919;
	volatility = 0.4251;
	flavor = Option::Call;
	discount = yield_curve->discount(T);
	strikepayoff = ext::make_shared<PlainVanillaPayoff>(flavor, strike);
	stddev = volatility * sqrt(T);
	black = BlackCalculator(strikepayoff, spot, stddev, discount);
	cout << "Option Price: " << black.value() << endl;
	cout << "Delta: " << black.delta(spot) << endl;
	cout << "Gamma: " << black.gamma(spot) << endl;
	cout << "Theta: " << black.theta(spot, t) << endl;
	cout << "Vega: " << black.vega(t) << endl;
	cout << "Rho: " << black.rho(t) << endl;

}

Real greek(const VanillaOption &option ,const ext::shared_ptr<SimpleQuote> &quote, const Real &dx) { 
	Real x0 = quote->value();
	quote->setValue(x0 + dx);
	Real P_u = option.NPV();
	quote->setValue(x0 - dx);
	Real P_d = option.NPV();
	quote->setValue(x0);
	return (P_u - P_d) / (2 * dx);
}
void rho_for_the_black_process() { 
	Date today(24, December, 2016);
	Settings::instance().evaluationDate() = today;
	auto u = ext::make_shared<SimpleQuote>(100.0);
	auto r = ext::make_shared<SimpleQuote>(0.01);
	auto sigma = ext::make_shared<SimpleQuote>(0.20);

	auto risk_free_curve = ext::make_shared<FlatForward>(today, Handle<Quote>(r), Actual360());
	auto dividend_ts_curve = ext::make_shared<FlatForward>(today, Handle<Quote>(r), Actual360());
	auto volatility = ext::make_shared<BlackConstantVol>(today, TARGET(), Handle<Quote>(sigma), Actual360());

	auto process_1 = ext::make_shared<BlackScholesProcess>(Handle<Quote>(u), Handle<YieldTermStructure>(risk_free_curve), Handle<BlackVolTermStructure>(volatility));
	cout << process_1->dividendYield()->zeroRate(1.0, Continuous).rate() << endl;

	auto process_2 = ext::make_shared<BlackProcess>(Handle<Quote>(u), Handle<YieldTermStructure>(risk_free_curve), Handle<BlackVolTermStructure>(volatility));
	cout << process_2->riskFreeRate()->zeroRate(1.0, Continuous).rate() << endl;
	cout << process_2->dividendYield()->zeroRate(1.0, Continuous).rate() << endl;


	EuropeanOption option_1(ext::make_shared<PlainVanillaPayoff>(Option::Call, 100.0), ext::make_shared<EuropeanExercise>(today + 100));
	option_1.setPricingEngine(ext::make_shared<AnalyticEuropeanEngine>(process_1));
	cout << option_1.NPV() << endl;

	EuropeanOption option_2(ext::make_shared<PlainVanillaPayoff>(Option::Call, 100.0), ext::make_shared<EuropeanExercise>(today + 100));
	option_2.setPricingEngine(ext::make_shared<AnalyticEuropeanEngine>(process_2));
	cout << option_2.NPV() << endl;

	cout << option_1.delta() << "  " << greek(option_1, u, 0.01) << endl;
	cout << option_2.delta() << "  " << greek(option_2, u, 0.01) << endl;

	cout << option_1.vega() << "  " << greek(option_1, sigma, 0.01) << endl;
	cout << option_2.vega() << "  " << greek(option_2, sigma, 0.01) << endl;

	cout << option_1.rho() << "  " << greek(option_1, r, 0.01) << endl;
	cout << option_2.rho() << "  " << greek(option_2, r, 0.01) << endl;

	cout << option_2.rho()+option_2.dividendRho() << "  " << greek(option_2, r, 0.001) << endl;

	auto process_3 = ext::make_shared<BlackScholesMertonProcess>(Handle<Quote>(u), Handle<YieldTermStructure>(risk_free_curve), Handle<YieldTermStructure>(dividend_ts_curve),
										Handle<BlackVolTermStructure>(volatility));
	EuropeanOption option_3(ext::make_shared<PlainVanillaPayoff>(Option::Call, 100.0), ext::make_shared<EuropeanExercise>(today+100));
	option_3.setPricingEngine(ext::make_shared<AnalyticEuropeanEngine>(process_3));

	cout << option_3.delta() << " " << greek(option_3, u, 0.01) << endl;
	cout << option_3.rho() << " " << greek(option_3, r, 0.001) << "  "<< option_3.rho() + option_3.dividendRho() << endl;
}