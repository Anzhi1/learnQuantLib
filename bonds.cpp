#include <iostream>
#include <ql/quantlib.hpp>
#include "bonds.h"

using namespace std;
using namespace QuantLib;

void fixed_rate_bond() {
	
	Real npv = 3 / (pow(1+0.005,0.5)) + (100+3)/(1+0.007);
	cout << npv << endl;

	Date calc_date(15, January, 2015);
	Settings::instance().evaluationDate() = calc_date;
	vector<Date> spot_dates = { Date(15, January, 2015), Date(15, July, 2015), Date(15, January, 2016) };
	vector<Rate> spot_rates = { 0.0, 0.005, 0.007 };
	DayCounter day_count = Thirty360(Thirty360::BondBasis);
	Calendar calendar = UnitedStates(UnitedStates::Settlement);
	auto interpolation = Linear();
	Compounding compounding = Compounded;
	Frequency compounding_frequency = Annual;
	auto spot_curve = ext::make_shared<ZeroCurve>(spot_dates, spot_rates, day_count, calendar, interpolation, compounding, compounding_frequency);
	Handle<YieldTermStructure> spot_curve_handle(spot_curve); 

	Date issue_Date(15, January, 2015);
	Date maturity_date(15, January, 2016);
	Period tenor(Semiannual);
	BusinessDayConvention business_convention = Unadjusted;
	DateGeneration::Rule rule = DateGeneration::Backward;
	bool month_end = false;
	Schedule schedule(issue_Date, maturity_date, tenor, calendar, business_convention,business_convention, rule, month_end);
	auto dates = schedule.dates();
	for (auto& date:dates) {
		cout << date << " ";
	}

	Real coupon_rate = .06;
	vector<Real> coupons = { coupon_rate };
	Natural settlement_days = 0;
	Real face_value = 100;
	FixedRateBond fixed_rate_bond(settlement_days, face_value, schedule, coupons, day_count);
	auto bond_engine = ext::make_shared<DiscountingBondEngine>(spot_curve_handle);
	fixed_rate_bond.setPricingEngine(bond_engine);
	cout << endl;
	cout << fixed_rate_bond.NPV() << endl;
	cout << fixed_rate_bond.cleanPrice() << endl; 
	cout << fixed_rate_bond.accruedAmount() << endl;
	cout << fixed_rate_bond.dirtyPrice() << endl; 
//	cout << BondFunctions::yield(fixed_rate_bond, face_value, day_count, compounding, compounding_frequency) << endl;
	cout << fixed_rate_bond.yield(day_count, compounding, compounding_frequency);
	
}

CallableFixedRateBond value_bond(const Handle<YieldTermStructure> &ts_handle,const Real &a, const Real &s, const Natural &grid_points, CallableFixedRateBond &bond) {
	auto model = ext::make_shared<HullWhite>(ts_handle);
	auto engine = ext::make_shared<TreeCallableFixedRateBondEngine>(model,grid_points);
	bond.setPricingEngine(engine);
	return bond; 
}
void callable_bond() {
	Date calc_date(16, August, 2016);
	Settings::instance().evaluationDate() = calc_date;
	DayCounter day_count = ActualActual(ActualActual::Bond);
	Real rate = 0.035;
	auto ts = ext::make_shared<FlatForward>(calc_date, rate, day_count, Compounded, Semiannual);
	Handle<YieldTermStructure> ts_handle(ts);

	CallabilitySchedule callability_schedule;
	Real call_price = 100.0;
	Date call_date(15, September, 2016);
	Calendar null_calendar = NullCalendar();

	for (int i = 0; i < 24; i++) {
		Bond::Price callability_price(call_price, Bond::Price::Clean);
		callability_schedule.push_back(ext::make_shared<Callability>(callability_price, Callability::Call,call_date));
		call_date = null_calendar.advance(call_date, 3, Months);
	}

	Date issue_date(16, September, 2014);
	Date maturity_date(15, September, 2022);
	Calendar calendar = UnitedStates(UnitedStates::GovernmentBond);
	Period tenor(Quarterly);
	BusinessDayConvention accrual_convention = Unadjusted;
	Schedule schedule(issue_date, maturity_date, tenor, calendar, accrual_convention, accrual_convention, DateGeneration::Backward, false);
	Natural settlement_days = 3;
	Real face_amount = 100;
	DayCounter accrual_daycount = ActualActual(ActualActual::Bond);
	Real coupon = 0.025;
	CallableFixedRateBond bond(settlement_days, face_amount, schedule, {coupon}, accrual_daycount, Following, face_amount, issue_date, callability_schedule);
	auto model = ext::make_shared<HullWhite>(ts_handle, 0.03, 0.12);
	auto engine = ext::make_shared<TreeCallableFixedRateBondEngine>(model,40);
	bond.setPricingEngine(engine);
//	value_bond(ts_handle, 0.03, 0.12, 40, bond);
	cout << "Bond price: " << bond.NPV() << endl;


    auto cashflows = bond.cashflows();
	for (int i = 0; i < cashflows.size(); i++) {
		cout << cashflows[i]->date() << "  " << cashflows[i]->amount() << endl;
	}
/*
	cout << "--------callability schedule----------" << endl;
	auto callability = bond.callability();
	for (int i = 0; i < callability.size(); i++) {
		cout << callability[i]->date() << "  " << callability[i]->price().amount() << endl;
	}
*/


}

void duration_of_floating_rate_bonds() {
	Date today(8, October, 2014);
	Settings::instance().evaluationDate() = today;

	RelinkableHandle<YieldTermStructure> forecast_curve;
	forecast_curve.linkTo(ext::make_shared<FlatForward>(today, 0.002, Actual360(),Compounded, Semiannual));
	auto index = ext::make_shared<Euribor6M>(forecast_curve);
	index->addFixing(Date(6, August, 2014), 0.002);

	Date issueDate(8, August, 2014);
	Date maturityDate(8, August, 2019);
	Schedule schedule(issueDate, maturityDate, Period(Semiannual), TARGET(), Following, Following, DateGeneration::Backward, false);
	FloatingRateBond bond(3,100,schedule,index,Actual360());
	auto cashflows = bond.cashflows();
	auto y2 = ext::make_shared<SimpleQuote>(0.002);
	auto yield_curve = ext::make_shared<FlatForward>(bond.settlementDate(), Handle<Quote>(y2), Actual360(), Compounded, Semiannual);
	vector<Date> dates;
	vector<Real> cfs;
	vector<Real> discounts;

	for (auto& cashflow : cashflows) {
		cout << cashflow->date() << "  " << cashflow->amount() << endl;
		dates.push_back(cashflow->date());
		cfs.push_back(cashflow->amount());
		discounts.push_back(yield_curve->discount(cashflow->date()));

	}
	
	InterestRate y(0.002, Actual360(), Compounded, Semiannual); 
	cout << BondFunctions::duration(bond, y, Duration::Modified) << endl;

	Real p = 0;
	for (int i = 0; i < cfs.size(); i++) {
		p += cfs[i] * discounts[i];
	}
	cout << "p: " << p << endl;
	bond.setPricingEngine(ext::make_shared<DiscountingBondEngine>(Handle<YieldTermStructure>(yield_curve)));
	cout << bond.dirtyPrice() << endl;

	Real dy = 1e-5;
	y2->setValue(0.002 + dy);
	Real P_p = 0;
	for (auto cf : bond.cashflows()) {
		P_p += cf->amount() * yield_curve->discount(cf->date());
	}
	cout << "P_p: " << P_p << endl;

	y2->setValue(0.002 - dy);
	Real P_m = 0;
	for (auto cf : bond.cashflows()) {
		P_m += cf->amount() * yield_curve->discount(cf->date());
	}
	cout << "P_m: " << P_m << endl;
	cout << -(1 / p) * (P_p - P_m) / (2 * dy) << endl;

	//the problem above is that it doesn't use the yield curve for forecasting, so it's not really considering the bond as a floating-rate bond,
	//it's using it as a fixed-rate bond, whose coupon rates happen to equal the current forecasts for the Euribor 6M fixings. This is clear if we
	//look at the coupon amounts and discounts we stored during the calculations
	//the discount factors changed when the yield was modified, but the coupon amounts stayed the same.
	cout << "--------------link the forecast_curve to the yield_curve we modified---------------" << endl;
	forecast_curve.linkTo(yield_curve);
	y2->setValue(0.002 + dy);
	P_p = bond.dirtyPrice();

	cout << "P_p: " << P_p << endl;

	y2->setValue(0.002 - dy);
	P_m = bond.dirtyPrice();
	y2->setValue(0.002);

	cout << "P_m: " << P_m << endl;
	cout << -(1 / p) * (P_p - P_m) / (2 * dy) << endl;

	auto discount_curve = ext::make_shared<ZeroSpreadedTermStructure>(forecast_curve, Handle<Quote>(ext::make_shared<SimpleQuote>(0.001)));
	bond.setPricingEngine(ext::make_shared<DiscountingBondEngine>(Handle<YieldTermStructure>(discount_curve)));

	p = bond.dirtyPrice();
	cout << "after we changed the discount_curve, p:" << p << endl;

	y2->setValue(0.002 + dy);
	P_p = bond.dirtyPrice();
	cout << "P_p: " << P_p << endl;

	y2->setValue(0.002 - dy);
	P_m = bond.dirtyPrice();
	cout << "P_m: " << P_m << endl;
	y2->setValue(0.002);
	cout << -(1 / p) * (P_p - P_m) / (2 * dy) << endl;


}

FixedRateBond create_tsy_security(const Date &bond_issue_date, const Date &bond_maturity_date, const Real &coupon_rate, Period coupon_frequency = Period(6, Months),
	DayCounter day_count = ActualActual(ActualActual::Bond), Calendar calendar = UnitedStates(UnitedStates::GovernmentBond)) {
	
	Real face_value = 100;
	Natural settlement_days = 0;

	Schedule schedule(bond_issue_date, bond_maturity_date, coupon_frequency, calendar, ModifiedFollowing, ModifiedFollowing, DateGeneration::Backward, false);
	FixedRateBond security(settlement_days, face_value, schedule, vector<Real>{coupon_rate}, day_count);
	return security;

}

void treasury_futures_contact() {
	Date calc_date(30, November, 2015);
	Settings::instance().evaluationDate() = calc_date;
	DayCounter day_count = ActualActual(ActualActual::Bond);
	Calendar calendar = UnitedStates(UnitedStates::Settlement);
	BusinessDayConvention business_convention = Following;
	bool end_of_month = false;
	Natural settlement_days = 0;
	Real face_amount = 100;
	Period coupon_frequency(Semiannual);
	
	vector<Real> prices = { 99.9935, 99.9576, 99.8119, 99.5472, 99.8867, 100.0664, 99.8711, 100.0547, 100.3047, 100.2266 };
	vector<Real> coupon_rates = { 0.0000, 0.0000, 0.0000, 0.0000, 0.00875, 0.0125, 0.01625, 0.02, 0.0225, 0.03 };

	vector<Date> maturity_dates = { Date(24,December,2015), Date(25,February,2016), Date(26,May,2016), Date(10,November,2016), Date(30,November,2017),
		Date(15,November,2018), Date(30,November,2020), Date(30,November,2022), Date(15,November,2025), Date(15,November,2045) };

	vector<ext::shared_ptr<RateHelper>> bond_helpers;
	
	for (int i = 0; i < prices.size(); i++) { 
		Schedule schedule(calc_date, maturity_dates[i], coupon_frequency, calendar, business_convention, business_convention,
			DateGeneration::Backward, false);

		auto helper = ext::make_shared<FixedRateBondHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(prices[i])), settlement_days, face_amount, schedule,
			vector<Real>{coupon_rates[i]}, day_count, business_convention);
		bond_helpers.push_back(helper);
	}

	auto yield_curve = ext::make_shared<PiecewiseYieldCurve<ZeroYield, Cubic>>(calc_date, bond_helpers, day_count);
	Handle<YieldTermStructure> yield_curve_handle(yield_curve);
	vector<Real> discount_factors;
	for (auto& date : maturity_dates) {
		Real discount = yield_curve->discount(date);
		discount_factors.push_back(discount); 
		cout << date << "  " << discount << endl;
	}


	Date bond_issue_date = calc_date;
	Date delivery_date(1, December, 2015);
	Date bond_maturity_date = bond_issue_date + Period(10, Years);
	day_count = ActualActual(ActualActual::Bond);
	Real coupon_rate = 6 / 100.0;

	FixedRateBond deliverable = create_tsy_security(bond_issue_date,bond_maturity_date, coupon_rate, coupon_frequency, day_count, calendar);
	auto bond_engine = ext::make_shared<DiscountingBondEngine>(yield_curve_handle);
	deliverable.setPricingEngine(bond_engine);

	//lets calculate the Z-Spread for this deliverable .The Z-Spread is the static spread added to the yield curve to match the price of the security
	//This spread is a measure of the risk associated with the security .For a treasury security , you would expect this to be zero
	 
	Real futures_price = 127.0625;
	Real clean_price = futures_price * yield_curve->discount(delivery_date);
	Real zspread = BondFunctions::zSpread(deliverable, clean_price, yield_curve, day_count, Compounded, Semiannual, calc_date)*10000;
	cout << "z-spread: " << zspread << "bp" << endl;


	day_count = ActualActual(ActualActual::Bond);

	vector<double> coupons = { 1.625, 1.625, 1.75, 1.75, 1.875, 1.875, 2.0, 2.0, 2.0, 2.0, 2.125, 2.125, 2.25, 2.25, 2.375, 2.5, 2.5, 2.75, 2.75 };
	maturity_dates = {
		Date(15, August, 2022), Date(15, November, 2022), Date(30, September, 2022), Date(15, May, 2023), Date(31, August, 2022), 
		Date(31, October, 2022), Date(31, July, 2022), Date(15, February, 2023), Date(15, February, 2025), Date(15, August, 2025), 
		Date(30, June, 2022), Date(15, May, 2025), Date(15, November, 2024), Date(15, November, 2025), Date(15, August, 2024), 
		Date(15, August, 2023), Date(15, May, 2024), Date(15, November, 2023), Date(15, February, 2024) };

	
	prices = {
		97.921875, 97.671875, 98.546875, 97.984375, 99.375, 99.296875, 100.265625, 100.0625, 98.296875, 98.09375,
		101.0625, 99.25, 100.546875, 100.375, 101.671875, 103.25, 102.796875, 105.0625, 104.875
	};
	vector<FixedRateBond> securities;
	vector<Real> securities_prices;
	Real min_basis = 100, min_basis_index = -1;
	for (int i = 0; i < coupons.size(); i++) {
		Date issue_date = maturity_dates[i] - Period(10, Years);
		FixedRateBond s = create_tsy_security(issue_date, maturity_dates[i], coupons[i] / 100.0);
		auto bond_engine = ext::make_shared<DiscountingBondEngine>(yield_curve_handle);
		s.setPricingEngine(bond_engine);
		Real cf = BondFunctions::cleanPrice(s, 0.06, day_count, Compounded, Semiannual, calc_date)/100.0;
		Real adjusted_futures_price = futures_price * cf;
		Real basis = prices[i] - adjusted_futures_price;
		if (basis < min_basis) {
			min_basis = basis;
			min_basis_index = i;
		}
		securities.push_back(s);
		securities_prices.push_back(cf);
	}

	Real ctd_coupon = coupons[min_basis_index];
	Date ctd_maturity = maturity_dates[min_basis_index];
	Real ctd_price = prices[min_basis_index];
	Real ctd_cf = securities_prices[min_basis_index];
	FixedRateBond ctd_bond = securities[min_basis_index];
	cout << "Minimum Basis: " << min_basis << endl;
	cout << "Conversion Factor: " << ctd_cf << endl;
	cout << "Coupon " << ctd_coupon << endl;
	cout << "Maturity " << ctd_maturity << endl;
	cout << "Price " << ctd_price << endl;
	

	Date futures_maturity_date(21, December, 2015);
	BondForward futures(calc_date, futures_maturity_date, Position::Long, 0.0, settlement_days, day_count,
								 calendar, business_convention, ext::make_shared<FixedRateBond>(ctd_bond), yield_curve_handle, yield_curve_handle);
	Real model_futures_price = futures.cleanForwardPrice()/ctd_cf;
	auto implied_yield = futures.impliedYield(ctd_price/ctd_cf, futures_price, calc_date, Compounded, day_count).rate();
	Real z_spread = BondFunctions::zSpread(ctd_bond, ctd_price, yield_curve, day_count, Compounded, Semiannual, calc_date);
	auto ytm = BondFunctions::yield(ctd_bond, ctd_price, day_count, Compounded, Semiannual, calc_date);

	cout << "Model Futures Price " << model_futures_price << endl; 
	cout << "Market Futures Price " << futures_price << endl;
	cout << "Model Adjustment " << model_futures_price - futures_price << endl;
	cout << "Implied Yield " << implied_yield * 100 << endl;
	cout << "Forward Z-spread " << z_spread * 10000 << endl;
	cout << "Forward YTM " << ytm * 100 << endl;

}
void mischievous_pricing_conventions() {
	Date today(5, January, 2010);
	Settings::instance().evaluationDate() = today;
	RelinkableHandle<YieldTermStructure> discounting_curve;
	RelinkableHandle<YieldTermStructure> forecasting_curve;
	
	auto index = ext::make_shared<USDLibor>(Period(3,Months), forecasting_curve);
	Natural settlement_days = 0;
	Calendar calendar = NullCalendar();
	Real face_amount = 100.0;
	Schedule schedule(today, today+Period(4, Years), Period(3, Months), calendar, Unadjusted, Unadjusted, DateGeneration::Forward, false);
	FloatingRateBond bond(settlement_days, face_amount, schedule, index, Thirty360(Thirty360::BondBasis), Unadjusted, 0);
	bond.setPricingEngine(ext::make_shared<DiscountingBondEngine>(discounting_curve));

	auto flat_rate = ext::make_shared<FlatForward>(today, 0.10, Thirty360(Thirty360::BondBasis), Compounded, Quarterly);
	forecasting_curve.linkTo(flat_rate);
	discounting_curve.linkTo(flat_rate);

	cout << bond.cleanPrice() << endl; 

	cout << flat_rate->dayCounter() << endl;
	cout << ext::dynamic_pointer_cast<Coupon>(bond.cashflows()[0])->dayCounter() << endl;
	cout << index->dayCounter() << endl;

	auto cashflows = bond.cashflows();
	for (int i = 0; i < cashflows.size()-1; i++) {
		auto coupon = ext::dynamic_pointer_cast<Coupon>(cashflows[i]);
		cout << coupon->date() << " " << coupon->rate() << " " << coupon->accrualPeriod() << endl;
	}

	bond = FloatingRateBond(settlement_days, face_amount, schedule, index, Actual360(), Unadjusted, 0);
	bond.setPricingEngine(ext::make_shared<DiscountingBondEngine>(discounting_curve));

	auto flat_rate_2 = ext::make_shared<FlatForward>(today, 0.10, Actual360(), Compounded, Quarterly);
	forecasting_curve.linkTo(flat_rate_2);
	discounting_curve.linkTo(flat_rate_2);
	cout << "--------a little adjustment about conventions--------" << endl;
	cout << bond.cleanPrice() << endl;

	//There's still a small discrepancy which is likely due to date adjustments in the underlying USD libor fixings.The coupon rates are
	//much better overall, so we seem to be on the right track 

	cashflows = bond.cashflows();
	for (int i = 0; i < cashflows.size() - 1; i++) {
		auto coupon = ext::dynamic_pointer_cast<Coupon>(cashflows[i]);
		cout << coupon->date() << " " << coupon->rate() << " " << coupon->accrualPeriod() << endl;
	}

	//To get a (theoretical) par bond , we can use a custom index whose conventions match exactly those of the bond we wanted to use
	cout << "-----a custom index whose conventions match exactly those of the bond-----" << endl;
	auto index2 = ext::make_shared<IborIndex>("Mock Libor", Period(3,Months), 0, USDCurrency(), NullCalendar(), Unadjusted,
											false, Thirty360(Thirty360::BondBasis),forecasting_curve);
	bond = FloatingRateBond(settlement_days, face_amount, schedule, index2, Thirty360(Thirty360::BondBasis), Unadjusted, 0);
	bond.setPricingEngine(ext::make_shared<DiscountingBondEngine>(discounting_curve));
	forecasting_curve.linkTo(flat_rate);
	discounting_curve.linkTo(flat_rate); 
	cout << bond.cleanPrice() << endl;
	cashflows = bond.cashflows();
	for (int i = 0; i < cashflows.size() - 1; i++) {
		auto coupon = ext::dynamic_pointer_cast<Coupon>(cashflows[i]);
		cout << coupon->date() << " " << coupon->rate() << " " << coupon->accrualPeriod() << endl;
	}

}
struct Data {
	Real Amount;
    Time T;
    Real Discount;
	Real D_cumulative;
	Real A_discounted;
};


void more_mischievous_convention() {
	Date today(27, January, 2011);
	Settings::instance().evaluationDate() = today;

	Date issueDate(28, January, 2011);
	Date maturityDate(31, August, 2020);
	Schedule schedule(issueDate, maturityDate, Period(Semiannual), UnitedStates(UnitedStates::GovernmentBond), 
					Unadjusted, Unadjusted, DateGeneration::Backward, false);
	FixedRateBond bond(1, 100.0, schedule, vector<Real>{0.03625}, ActualActual(ActualActual::Bond), Unadjusted, 100.0);
	Real bond_yield = 0.034921;
	Real P1 = bond.dirtyPrice(bond_yield, bond.dayCounter(), Compounded, Semiannual);
	auto flat_curve = ext::make_shared<FlatForward>(bond.settlementDate(), bond_yield, ActualActual(ActualActual::Bond), Compounded, Semiannual); 
	auto engine = ext::make_shared<DiscountingBondEngine>(Handle<YieldTermStructure>(flat_curve));
	bond.setPricingEngine(engine);
	Real P2 = bond.dirtyPrice();
	cout << P1 << "  " << P2 << endl;

	DayCounter dayCounter = ActualActual(ActualActual::Bond);
	Time T = dayCounter.yearFraction(Date(28, January, 2011), Date(28, February, 2011), Date(28, August, 2010), Date(28, February, 2011));
	cout << T << endl;

	cout << dayCounter.yearFraction(Date(28, January, 2011), Date(28, February, 2011), Date(28, February, 2010), Date(28, February, 2011)) <<endl;

	InterestRate y(bond_yield, dayCounter, Compounded, Semiannual);
	cout << y.discountFactor(T) << endl;

	auto cashflows = bond.cashflows();
	vector<Data> data;
	Real D = 0;
	Real D_cumulative = 0 ;
	for (int i = 0; i < cashflows.size() - 1; i++) {
		auto c = ext::dynamic_pointer_cast<Coupon>(cashflows[i]);
		Real A = c->amount();
		Time T = c->accrualPeriod();
		D = y.discountFactor(T);
		D_cumulative = (i == 0 ? D : D * (*(data.end() - 1)).D_cumulative);
		Real A_discounted = A * D_cumulative; 
		data.push_back(Data{ A,T,D,D_cumulative,A_discounted });
	}
	data.push_back(Data{100,0,0,D_cumulative,100*D_cumulative});

	Real sum = 0;
	for (int i = 0; i < data.size(); i++) {
		cout << data[i].Amount << " " << data[i].T << " " << data[i].Discount << " " << data[i].D_cumulative << " " << data[i].A_discounted << endl;
		sum += data[i].A_discounted;
	}
	cout << sum << endl;
	
	cout << "----------curve-based calculation----------" << endl;

	vector<Real> amount;
	vector<Discount> discounts;
	vector<Real> a_discounts;
	sum = 0;
	for (int i = 0; i < cashflows.size(); i++) {
		Real A = cashflows[i]->amount();
		Real D = flat_curve->discount(cashflows[i]->date());
		Real A_discounted = A * D;
		cout << A << " " << D << " " << A_discounted << endl;
		sum += A_discounted;
	}
	cout << sum << endl;

	cout << bond.dirtyPrice() << endl;
	
	cout << flat_curve->discount(Date(28, February, 2011)) << endl;

	cout << dayCounter.yearFraction(Date(28, January, 2011), Date(28, February, 2011)) << endl;

	T = dayCounter.yearFraction(Date(28, January, 2011), Date(28, February, 2011), Date(28, August, 2010), Date(28, February, 2011));
	cout << T << endl;

	Real P_y = bond.dirtyPrice(bond_yield, bond.dayCounter(), Compounded, Semiannual);
	Real D_y = y.discountFactor(T);

	Real P_c = bond.dirtyPrice();
	Real D_c = flat_curve->discount(Date(28, February, 2011));

	cout << P_y << endl;
	cout << P_c * (D_y / D_c) <<endl;
}