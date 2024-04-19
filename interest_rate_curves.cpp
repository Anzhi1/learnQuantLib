#include "interest_rate_curves.h"
#include <ql/quantlib.hpp>
#include <iostream>
#include <cmath>

using namespace std;
using namespace QuantLib;

void term_structure_reference_date() {
	Settings::instance().evaluationDate() = Date(3, October, 2014);

	vector<ext::shared_ptr<RateHelper>> helpers;
	vector<int> tenors = {2, 3, 5, 10, 15};
	vector<Real> rates = { 0.201, 0.258, 0.464, 1.151, 1.588 };
	const auto index = ext::make_shared<Euribor6M>();
	for (int i = 0; i < tenors.size(); i++) {
		ext::shared_ptr<SimpleQuote> quote = ext::make_shared<SimpleQuote>(rates[i]/100.0);
		auto helper = ext::make_shared<SwapRateHelper>(Handle<Quote>(quote), Period(tenors[i], Years),
			TARGET(), Annual, Unadjusted, Thirty360(Thirty360::European), index);
		helpers.push_back(helper);
	}
	auto curve1 = ext::make_shared<PiecewiseYieldCurve<ForwardRate, BackwardFlat>>(0, TARGET(), helpers, Actual360());
	
	auto nodes = curve1->nodes();
	vector<Date> dates = {};
	rates.clear();
	for (int i = 0; i < nodes.size(); i++) {
		cout << nodes[i].first << " " << nodes[i].second << endl;
		dates.push_back(nodes[i].first);
		rates.push_back(nodes[i].second);
	}

	auto curve2 = ext::make_shared<InterpolatedForwardCurve<BackwardFlat>>(dates,rates,Actual360());
	
	cout << curve1->referenceDate() << " to " << curve1->maxDate() << endl;
	cout << curve2->referenceDate() << " to " << curve2->maxDate() << endl;

	cout << curve1->zeroRate(5.0, Continuous) << endl;
	cout << curve2->zeroRate(5.0, Continuous) << endl;

	cout << curve1->zeroRate(Date(7, September, 2019), Actual360(), Continuous) << endl;
	cout << curve2->zeroRate(Date(7, September, 2019), Actual360(), Continuous) << endl;

	//moving the evaluation date
	//to recap: we built the first curve specifying its reference date relative to the evaluation date, and the second curve specifying 
	//its reference date explicitly . Now what happens if we change the evaluation date

	Settings::instance().evaluationDate() = Date(19, September, 2014);

	//curve1 changes accordingly while curve2 doesn't
	cout << "after we changed the evaluation date" << endl;
	cout << curve1->referenceDate() << " to " << curve1->maxDate() << endl;
	cout << curve2->referenceDate() << " to " << curve2->maxDate() << endl;

	cout << curve1->zeroRate(5.0, Continuous) << endl;
	cout << curve2->zeroRate(5.0, Continuous) << endl;

	cout << curve1->zeroRate(Date(7, September, 2019), Actual360(), Continuous) << endl;
	cout << curve2->zeroRate(Date(7, September, 2019), Actual360(), Continuous) << endl;
}


void eonia_curve() {
	
	Date today(11, December, 2012);
	Settings::instance().evaluationDate() = today;
	vector<ext::shared_ptr<RateHelper>> helpers;
	Real rate = 0.04;
	for (int i = 0; i < 3; i++) {
		ext::shared_ptr<SimpleQuote> quote = ext::make_shared<SimpleQuote>(rate / 100);
		auto helper = ext::make_shared<DepositRateHelper>(Handle<Quote>(quote), Period(1,Days), i,TARGET(), Following, false, Actual360());
		helpers.push_back(helper);
	}

	auto eonia = ext::make_shared<Eonia>();
	vector<Rate> rates = { 0.070, 0.069, 0.078, 0.074 };
	vector<Period> tenors = { Period(1,Weeks), Period(2,Weeks), Period(3, Weeks),Period(1,Months) };
	for (int i = 0; i < rates.size(); i++) {
		auto quote = ext::make_shared<SimpleQuote>(rates[i] / 100);
		auto helper = ext::make_shared<OISRateHelper>(2, Period(tenors[i]),Handle<Quote>(quote),eonia);
		helpers.push_back(helper);
	}

	vector<Rate> datedOISrates = {0.046, 0.016, -0.007, -0.013, -0.014};
	vector<Date> startDates = { Date(16, January, 2013), Date(13, February, 2013), Date(13, March, 2013), Date(10, April, 2013), Date(8, May, 2013) };
	vector<Date> endDates =   {Date(13, February,2013),  Date(13, March, 2013),    Date(10, April, 2013), Date(8, May, 2013),    Date(12, June, 2013)};

	for (int i = 0; i < datedOISrates.size(); i++) {
		auto quote = ext::make_shared<SimpleQuote>(datedOISrates[i] / 100);
		auto helper = ext::make_shared<DatedOISRateHelper>(startDates[i], endDates[i], Handle<Quote>(quote), eonia);
		helpers.push_back(helper);
	}

	//finally we add OIS quotes up to 30 years
	vector<Rate> oisRates = {0.002, 0.008, 0.021, 0.036, 0.127, 0.274, 0.456, 0.647, 0.827, 0.996, 1.147, 1.280, 1.404, 1.516, 1.764, 1.939, 2.003, 2.038};
	vector<Period> oisTenor = { Period(15, Months), Period(18, Months), Period(21, Months), Period(2, Years), Period(3, Years), Period(4, Years), Period(5,Years),
	Period(6, Years), Period(7, Years), Period(8, Years), Period(9, Years), Period(10, Years), Period(11,Years), Period(12, Years),
	Period(15, Years), Period(20, Years), Period(25, Years),Period(30, Years) };

	for (int i = 0; i < oisRates.size(); i++) {
		auto quote = ext::make_shared<SimpleQuote>(oisRates[i] / 100);
		auto helper = ext::make_shared<OISRateHelper>(2, oisTenor[i], Handle<Quote>(quote), eonia);
		helpers.push_back(helper);
	}
	
	//according to the QuantLib-SWIG interface file the PiecewiseLogCubic in python refers to PiecewiseYieldCurve<Discount,MonotonicLogCubic> 
	//instead of PiecewiseYieldCurve<Discount,LogCubic> 
	
	auto eonia_curve_c = ext::make_shared<PiecewiseYieldCurve<Discount, MonotonicLogCubic>>(0, TARGET(), helpers, Actual365Fixed());
	eonia_curve_c->enableExtrapolation();

	//turn-of-year jumps

	auto eonia_curve_ff = ext::make_shared<PiecewiseYieldCurve<ForwardRate, BackwardFlat>>(0, TARGET(), helpers, Actual365Fixed());
	eonia_curve_ff->enableExtrapolation();

	auto nodes = eonia_curve_ff->nodes();
	vector<Date> tempDates;
	vector<Rate> tempRates;
	for (int i = 0; i < nodes.size(); i++) {
		if (i < 9) {
			cout << nodes[i].first << " " << nodes[i].second << endl;
		}
		tempDates.push_back(nodes[i].first);
		tempRates.push_back(nodes[i].second);
	}

	//to create a curve that doesn't include the jump, we replace the relevant forward rate with a simple average of the ones that precede and follow
	tempRates[6] = (tempRates[5] + tempRates[7]) / 2.0;
	auto tempCurve = ext::make_shared<ForwardCurve>(tempDates,tempRates,eonia_curve_ff->dayCounter());

	Date d1 = Date(31, December, 2012) - Period(1, Weeks);
	Date d2 = Date(31, December, 2012) + Period(1, Weeks);

	Rate F = eonia_curve_ff->forwardRate(d1, d2, Actual360(), Simple).rate();
	Rate F_1 = tempCurve->forwardRate(d1, d2, Actual360(), Simple).rate();

	cout << F << "  " << F_1  << endl;


	auto t12 = eonia_curve_ff->dayCounter().yearFraction(d1, d2);
	auto t_j = eonia_curve_ff->dayCounter().yearFraction(Date(31, December, 2012), Date(2, January, 2013));

	auto J = (F - F_1) * t12 / t_j;
	cout << J << endl;

	auto B = 1.0 / (1.0 + J * t_j);
	
	
	ext::shared_ptr<SimpleQuote> quote = ext::make_shared<SimpleQuote>(B);
	vector<Handle<Quote>> jumps = { Handle<Quote>(quote) };
//	Handle<Quote> test(quote);
//	jumps2.push_back(test);
	const vector<Date> jumpDates = { Date(31, December, 2012) };
	Natural settlementDays = 0;
	auto eonia_curve = ext::make_shared<PiecewiseYieldCurve<Discount, MonotonicLogCubic>> (settlementDays, TARGET(), helpers,Actual365Fixed(), jumps,jumpDates);
	
	eonia_curve->enableExtrapolation();

}

void euribor_curve() {

	Date today(11, December, 2012);
	Settings::instance().evaluationDate() = today;

	//discounting curve
	auto eonia = ext::make_shared<Eonia>();
	vector<ext::shared_ptr<RateHelper>> helpers;
	Real rate = 0.04;
	for (int i = 0; i < 3; i++) {
		ext::shared_ptr<SimpleQuote> quote = ext::make_shared<SimpleQuote>(rate / 100);
		auto helper = ext::make_shared<DepositRateHelper>(Handle<Quote>(quote), Period(1, Days), i, TARGET(), Following, false, Actual360());
		helpers.push_back(helper);
	}

	vector<Rate> rates = { 0.070, 0.069, 0.078, 0.074 };
	vector<Period> tenors = { Period(1,Weeks), Period(2,Weeks), Period(3, Weeks),Period(1,Months) };
	for (int i = 0; i < rates.size(); i++) {
		auto quote = ext::make_shared<SimpleQuote>(rates[i] / 100);
		auto helper = ext::make_shared<OISRateHelper>(2, Period(tenors[i]), Handle<Quote>(quote), eonia);
		helpers.push_back(helper);
	}

	vector<Rate> datedOISrates = { 0.046, 0.016, -0.007, -0.013, -0.014 };
	vector<Date> startDates = { Date(16, January, 2013), Date(13, February, 2013), Date(13, March, 2013), Date(10, April, 2013), Date(8, May, 2013) };
	vector<Date> endDates = { Date(13, February,2013),  Date(13, March, 2013),    Date(10, April, 2013), Date(8, May, 2013),    Date(12, June, 2013) };

	for (int i = 0; i < datedOISrates.size(); i++) {
		auto quote = ext::make_shared<SimpleQuote>(datedOISrates[i] / 100);
		auto helper = ext::make_shared<DatedOISRateHelper>(startDates[i], endDates[i], Handle<Quote>(quote), eonia);
		helpers.push_back(helper);
	}
	vector<Rate> oisRates = { 0.002, 0.008, 0.021, 0.036, 0.127, 0.274, 0.456, 0.647, 0.827, 0.996, 1.147, 1.280, 1.404, 1.516, 1.764, 1.939, 2.003, 2.038 };
	vector<Period> oisTenor = { Period(15, Months), Period(18, Months), Period(21, Months), Period(2, Years), Period(3, Years), Period(4, Years), Period(5,Years),
	Period(6, Years), Period(7, Years), Period(8, Years), Period(9, Years), Period(10, Years), Period(11,Years), Period(12, Years),
	Period(15, Years), Period(20, Years), Period(25, Years),Period(30, Years) };

	for (int i = 0; i < oisRates.size(); i++) {
		auto quote = ext::make_shared<SimpleQuote>(oisRates[i] / 100);
		auto helper = ext::make_shared<OISRateHelper>(2, oisTenor[i], Handle<Quote>(quote), eonia);
		helpers.push_back(helper);
	}

	vector<Handle<Quote>> jumps = {};
	jumps.push_back(Handle<Quote>(ext::make_shared<SimpleQuote>(exp(-0.00102 * 2.0 / 360))));
	jumps.push_back(Handle<Quote>(ext::make_shared<SimpleQuote>(exp(-0.00086 * 2.0 / 360))));
	vector<Date> jumpDates = { Date(31, December, 2012), Date(31, December, 2013) };
	
	auto eonia_curve = ext::make_shared<PiecewiseYieldCurve<Discount, MonotonicLogCubic>>(2, TARGET(), helpers, Actual365Fixed(), jumps,jumpDates);
	eonia_curve->enableExtrapolation();

	// 6-months Euribor, which is somewhat simpler due to having a number of quoted rates directly available for bootstrapping
	helpers.clear();
	//ugly c++
	helpers.push_back(ext::make_shared<DepositRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(0.312 / 100)), Period(6, Months), 3, TARGET(), Following, false, Actual360()));

	auto euribor6m = ext::make_shared<Euribor6M>();
	vector<Rate> fraRates = {0.293, 0.272, 0.260, 0.256, 0.252, 0.248, 0.254, 0.261, 0.267, 0.279, 0.291, 0.303, 0.318, 0.335, 0.352, 0.371, 0.389, 0.409};
	vector<int> start = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 };

	for (int i = 0; i < fraRates.size(); i++) {
		auto quote = ext::make_shared<SimpleQuote>(fraRates[i] / 100);
		auto helper = ext::make_shared<FraRateHelper>(Handle<Quote>(quote), start[i], euribor6m);
		helpers.push_back(helper);
	}

	RelinkableHandle<YieldTermStructure> discount_curve;
	discount_curve.linkTo(eonia_curve);

	vector<Rate> swapRates = { 0.424, 0.576, 0.762, 0.954, 1.135, 1.303, 1.452, 1.584, 1.809, 2.037, 2.187, 2.234, 2.256, 2.295, 2.348, 2.421, 2.463 };
	vector<int> tenor = {3,4,5,6,7,8,9,10,12,15,20,25,30,35,40,50,60 };

	for (int i = 0; i < swapRates.size(); i++) {
		auto quote = ext::make_shared<SimpleQuote>(swapRates[i]/100);
		auto helper = ext::make_shared<SwapRateHelper>(Handle<Quote>(quote), Period(tenor[i], Years), TARGET(), Annual,
			Unadjusted, Thirty360(Thirty360::BondBasis), euribor6m, Handle<Quote>(), Period(0, Days), discount_curve);
		helpers.push_back(helper);
	}

	auto euribor6m_curve = ext::make_shared<PiecewiseYieldCurve<Discount, MonotonicLogCubic>>(2, TARGET(), helpers, Actual365Fixed());
	euribor6m_curve->enableExtrapolation();


	//Synthetic deposits

	ext::shared_ptr<PiecewiseYieldCurve<Discount, MonotonicLogCubic>> euribor6m_curve_0(euribor6m_curve);

	Date spot = euribor6m_curve->referenceDate();
	Date d = TARGET().advance(spot, 1, Days);

	auto F_x = euribor6m_curve_0->forwardRate(d, TARGET().advance(d, 6, Months), Actual360(), Simple).rate();
	auto F_on = eonia_curve->forwardRate(d, TARGET().advance(d, 6, Months), Actual360(), Simple).rate();
	DayCounter day_counter = euribor6m->dayCounter();
	auto T_x = day_counter.yearFraction(d, TARGET().advance(d, 6, Months));
	Real alpha = F_x - F_on;
	cout << "alpha: " << alpha << endl;

	//from the basis , we can instantiate synthetic deposits for a number of maturities below 6 months;

	vector<ext::shared_ptr<RateHelper>> synth_helpers;

	vector<Period> periods = {Period(1, Days), Period(1, Weeks), Period(2, Weeks), Period(3, Weeks), Period(1, Months), Period(2, Months),
							  Period(3, Months), Period(4, Months), Period(5, Months) };
	for (int i = 0; i < periods.size(); i++) {
		Integer n = periods[i].length();
		TimeUnit units = periods[i].units();
		Time t = day_counter.yearFraction(spot, TARGET().advance(spot, n, units));
		F_on = eonia_curve->forwardRate(spot, TARGET().advance(spot, n, units), Actual360(), Simple).rate();

		Real F = F_on + alpha;
		cout << periods[i] << " " << F << endl;
		synth_helpers.push_back(ext::make_shared <DepositRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(F)),
			periods[i], 2, TARGET(), Following, false, Actual360()));
	}
	
	vector<ext::shared_ptr<RateHelper>> temphelpers(helpers);
	temphelpers.insert(temphelpers.end(), synth_helpers.begin(), synth_helpers.end());
    euribor6m_curve = ext::make_shared<PiecewiseYieldCurve<Discount, MonotonicLogCubic>>(2, TARGET(), temphelpers, Actual365Fixed());
	euribor6m_curve->enableExtrapolation();


	//if we wanted to determine more coefficients for the basis, we'd have to select more quotes and solve a linear system,
	//For instance , to determine both alpha and beta , we can use the TOM 6-months and the 1x7FRAS:
	Date startDate = TARGET().advance(spot, 1, Days);
	Date endDate = TARGET().advance(startDate, 6, Months);
	F_x = euribor6m_curve_0->forwardRate(startDate, endDate, Actual360(), Simple).rate();
	F_on = eonia_curve->forwardRate(startDate, endDate, Actual360(), Simple).rate();
	Real T_x0 = day_counter.yearFraction(startDate, endDate);
	Real Delta0 = F_x - F_on;

	startDate = TARGET().advance(spot, 1, Months);
	endDate = TARGET().advance(startDate, 6, Months);
	F_x = euribor6m_curve_0->forwardRate(startDate, endDate, Actual360(), Simple).rate();
	F_on = eonia_curve->forwardRate(startDate, endDate, Actual360(), Simple).rate();
	Real T_x1 = day_counter.yearFraction(startDate, endDate);
	Real Delta1 = F_x - F_on;
	Time t1 = day_counter.yearFraction(spot, startDate);
	Time t2 = day_counter.yearFraction(spot, endDate);


	Matrix L(2, 2);
	L[0][0] = T_x0;
	L[0][1] = 0.5 * std::pow(T_x0, 2);
	L[1][0] = T_x1;
	L[1][1] = 0.5 * (std::pow(t2, 2) - std::pow(t1, 2));
	
	Array b(2);
	b[0] = Delta0 * T_x0;
	b[1] = Delta1 * T_x1;

	//solve the linear equations
	Array solution = inverse(L) * b;
	alpha = solution[0];
	Real beta = solution[1];

	cout << "alpha: " << alpha << " beta: " << beta << endl;

	//again we can create the synthetic deposits
	synth_helpers.clear();
	for (int i = 0; i < periods.size(); i++) {
		Integer n = periods[i].length();
		TimeUnit units = periods[i].units();
		Time t = day_counter.yearFraction(spot, TARGET().advance(spot, n, units));
		F_on = eonia_curve->forwardRate(spot, TARGET().advance(spot, n, units), Actual360(), Simple).rate();
		Real F = F_on + alpha + 0.5 * beta * t;
		cout << periods[i] << " " << F << endl; 
		synth_helpers.push_back(ext::make_shared <DepositRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(F)), 
			periods[i], 2, TARGET(), Following, false, Actual360())); 
	} 
	temphelpers = helpers;
	temphelpers.insert(temphelpers.end(), synth_helpers.begin(), synth_helpers.end()); 
	euribor6m_curve = ext::make_shared<PiecewiseYieldCurve<Discount, MonotonicLogCubic>>(2, TARGET(), temphelpers, Actual365Fixed()); 
	euribor6m_curve->enableExtrapolation(); 

	//12-month euribor
	auto euribor12m = ext::make_shared<Euribor1Y>();
	helpers.clear();
	helpers.push_back(ext::make_shared<DepositRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(0.54 / 100)), Period(12, Months),
		2, TARGET(), Following, false, Actual360()));
	helpers.push_back(ext::make_shared<FraRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(0.5070 / 100)), 12, euribor12m));
	rates = { 0.424, 0.576, 0.762, 0.954, 1.135, 1.303, 1.452, 1.584, 1.809, 2.037, 2.187, 2.234, 2.256 };
	vector<Real> basis = { 0.179, 0.164, 0.151, 0.139, 0.130, 0.123, 0.118, 0.113, 0.106, 0.093, 0.080, 0.072, 0.066 };
	tenor = { 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30 };
	for (int i = 0; i < rates.size(); i++) {
		helpers.push_back(ext::make_shared<SwapRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>((rates[i] + basis[i]) / 100)),
			Period(tenor[i], Years), TARGET(), Annual, Unadjusted, Thirty360(Thirty360::BondBasis),
			euribor12m, Handle<Quote>(), Period(0, Days), discount_curve));
	}

	//again we'll be using synthetic helpers to improve the shape of the short end of the curve , The same procedure we used for the Euribor 6M curve lets us
	//create deposits with a number of maturities below 1 year. I'll skip the calculation and just create helpers with the resulting rates as reported by the
	//paper

	synth_helpers.clear();
	rates = { 0.6537, 0.6187, 0.5772, 0.5563 };
	tenor = { 1, 3, 6, 9 };
	for (int i = 0; i < rates.size(); i++) {
		synth_helpers.push_back(ext::make_shared<DepositRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(rates[i]/100)),
			Period(tenor[i],Years),2, TARGET(), Following, false, Actual360()));
	}
	rates = { 0.4974, 0.4783, 0.4822, 0.5481, 0.6025 };
	vector<int> months_to_start = { 3, 6, 9, 15, 18 };

	Real last_basis = 0.066;
	rates = { 2.295, 2.348, 2.421, 2.463 };
	tenor = { 35, 40, 50, 60 };
	for (int i = 0; i < rates.size(); i++) {
		synth_helpers.push_back(ext::make_shared<SwapRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>((rates[i]+last_basis)/100)),
			Period(tenor[i], Years),TARGET(), Annual, Unadjusted, Thirty360(Thirty360::BondBasis), euribor12m, Handle<Quote>(), Period(0, Days), discount_curve));
	}
	temphelpers = helpers;
	temphelpers.insert(temphelpers.end(),synth_helpers.begin(),synth_helpers.end());

	auto euribor12m_curve = ext::make_shared<PiecewiseYieldCurve<Discount, MonotonicLogCubic>>(2, TARGET(), temphelpers, Actual365Fixed()); 
	euribor12m_curve->allowsExtrapolation();

	auto euribor12m_curve_0 = ext::make_shared<PiecewiseYieldCurve<Discount, MonotonicLogCubic>>(2, TARGET() ,temphelpers, Actual365Fixed());
	euribor12m_curve_0->allowsExtrapolation();



	//3-months Euribor
	auto euribor3m = ext::make_shared<Euribor3M>();
	helpers.clear();
	helpers.push_back(ext::make_shared<DepositRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(0.179 / 100)), Period(3, Months), 3, TARGET(),
		Following, false, Actual360()));

	rates = { 0.1775, 0.1274, 0.1222, 0.1269, 0.1565, 0.1961, 0.2556, 0.3101 };
	vector<Date> dates = { Date(19, December, 2012), Date(20, March, 2013), Date(19, June, 2013), Date(18, September, 2013),
								  Date(18, December, 2013), Date(19, March, 2014), Date(18, June, 2014), Date(17, September, 2014) };

	for (int i = 0; i < rates.size(); i++) {
		auto helper = ext::make_shared<FuturesRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(100-rates[i])), dates[i], euribor3m, Handle<Quote>());
		helpers.push_back(helper);
	}

	rates = { 0.424, 0.576, 0.762, 0.954, 1.135, 1.303, 1.452, 1.584, 1.809, 2.037, 2.187, 2.234, 2.256, 2.348, 2.421 };
	basis = { 0.1395, 0.1390, 0.1395, 0.1375, 0.1350, 0.1320, 0.1285, 0.1250, 0.1170, 0.1045, 0.0885, 0.0780, 0.0700, 0.0600, 0.0540 };
	tenor = { 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30, 40, 50 };

	for (int i = 0; i < rates.size(); i++) {
		auto helper = ext::make_shared<SwapRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>((rates[i]-basis[i])/100)),Period(tenor[i], Years),
			TARGET(), Annual, Unadjusted, Thirty360(Thirty360::BondBasis), euribor3m, Handle<Quote>(), Period(0, Days),discount_curve);
		helpers.push_back(helper);
	}

	synth_helpers.clear();

	rates = { 0.1865, 0.1969, 0.1951, 0.1874 };
	periods = { Period(2, Weeks), Period(3, Weeks), Period(1, Months), Period(2, Months) };
	
	for (int i = 0; i < rates.size(); i++) {
		synth_helpers.push_back(ext::make_shared<DepositRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(rates[i]/100)), periods[i], 2, TARGET(), 
			Following, false, Actual360()));
	}

	rates = { 2.295, 2.463 };
	basis = { 0.0650, 0.0540 };
	tenor = { 35, 60 };
	for (int i = 0; i < rates.size(); i++) {
		synth_helpers.push_back(ext::make_shared<SwapRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>((rates[i] - basis[i])/100)),Period(tenor[i], Years),
			TARGET(),Annual, Unadjusted, Thirty360(Thirty360::BondBasis),euribor3m, Handle<Quote>(), Period(0, Days), discount_curve));
	}

	//turn of the year

	vector<Rate> futureQuotes = { 0.1775, 0.1274, 0.1222, 0.1269, 0.1565, 0.1961, 0.2556, 0.3101 };
	vector<Date> futureDates = {Date(19, December, 2012), Date(20, March, 2013), Date(19, June, 2013), Date(18, September, 2013), Date(18, December, 2013),
		Date(19, March , 2014), Date(18, June, 2014), Date(17, September, 2014)};

	spot = euribor6m_curve->referenceDate();
	day_counter = euribor3m->dayCounter();
	vector<Time> times;
	Time t4= 0;
	for (int i = 0; i < futureDates.size();i++) {
		if (i != 4) {
			times.push_back(day_counter.yearFraction(spot, futureDates[i]));
		}
		else {
			t4 = day_counter.yearFraction(spot, futureDates[i]);
		}
		
	}
	vector<Rate> tempRates = { 0.1775, 0.1274, 0.1222, 0.1269, 0.1961, 0.2556, 0.3101 };
	MonotonicCubicNaturalSpline f(times.begin(),times.end(), tempRates.begin());
	cout << futureQuotes[4] << " "  << " " << f(t4) << endl;


	Real J = (futureQuotes[4] - f(t4)) / 100;
	Real tau = day_counter.yearFraction(Date(18, December, 2013), Date(18, March, 2014));
	cout << J <<"  " << tau << endl;

	jumps = { Handle<Quote>(ext::make_shared<SimpleQuote>(exp(-J * tau))) };
	jumpDates = { Date(31, December, 2012) };

	temphelpers = helpers;
	temphelpers.insert(temphelpers.end(),synth_helpers.begin(),synth_helpers.end());
	auto euribor3m_curve = ext::make_shared<PiecewiseYieldCurve<Discount, MonotonicLogCubic>>(2, TARGET(), temphelpers, Actual365Fixed(), jumps, jumpDates);
	euribor3m_curve->enableExtrapolation();



	// 1-month Euribor
	auto euribor1m = ext::make_shared<Euribor1M>();
	helpers.clear();
	synth_helpers.clear();
	helpers.push_back(ext::make_shared<DepositRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(0.110/100)),Period(1, Months), 2, TARGET(), Following,
		false, Actual360()));
	
	rates = { 0.106, 0.096, 0.085, 0.079, 0.075, 0.071, 0.069, 0.066, 0.065, 0.064, 0.063 };
	tenor = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
	for (int i = 0; i < rates.size(); i++) {
		helpers.push_back(ext::make_shared<SwapRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(rates[i] / 100)), Period(tenor[i], Months), TARGET(), Monthly,
			Unadjusted, Thirty360(Thirty360::BondBasis), euribor1m, Handle<Quote>(), Period(0, Days), discount_curve));
	}
	// for longer maturities , we can combine the swaps against 6-months euribor with the 1-month vs 6-months basis swaps 

	rates = { 0.324, 0.424, 0.576, 0.762, 0.954, 1.135, 1.303, 1.452, 1.584, 1.703, 1.809, 2.037, 2.187, 2.234, 2.256 };
	basis = { 0.226, 0.238, 0.246, 0.250, 0.250, 0.248, 0.245, 0.241, 0.237, 0.233, 0.228, 0.211, 0.189, 0.175, 0.163 };
	tenor = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 20, 25, 30 };
	for (int i = 0; i < rates.size(); i++) {
		helpers.push_back(ext::make_shared<SwapRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>((rates[i]-basis[i])/100)), Period(tenor[i], Years), 
			TARGET(), Annual, Unadjusted, Thirty360(Thirty360::BondBasis), euribor1m, Handle<Quote>(), Period(0, Days), discount_curve));
	}
	//As before, we can use synthetic deposits for maturities below the 1-month tenor
	rates = { 0.0661, 0.098, 0.0993, 0.1105 };
	periods = { Period(1, Days), Period(1, Weeks), Period(2, Weeks), Period(3, Weeks) };

	for (int i = 0; i < rates.size(); i++) {
		synth_helpers.push_back(ext::make_shared<DepositRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(rates[i]/100)), periods[i],
			2, TARGET(), Following, false, Actual360()));
	}

	last_basis = 0.163;
	rates = { 2.295, 2.348, 2.421, 2.463 };
	tenor = { 35, 40, 50, 60 };
	for (int i = 0; i < rates.size(); i++) {
		synth_helpers.push_back(ext::make_shared<SwapRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>((rates[i]-last_basis)/100)), 
			Period(tenor[i],Years),TARGET(), Annual, Unadjusted, Thirty360(Thirty360::BondBasis), euribor1m, Handle<Quote>(), Period(0,Days), discount_curve));
	}

	J = 0.0016;
	auto t_j = euribor1m->dayCounter().yearFraction(Date(31, December, 2012), Date(2,January, 2013));
	Real B = 1.0 / (1.0 + J * t_j);
	jumps = { Handle<Quote>(ext::make_shared<SimpleQuote>(B)) };
	jumpDates = {Date(31, December, 2013)};

	temphelpers = helpers;
	temphelpers.insert(temphelpers.end(), synth_helpers.begin(), synth_helpers.end());
	auto euribor1m_curve = ext::make_shared<PiecewiseYieldCurve<Discount, MonotonicLogCubic>>(2, TARGET(), temphelpers, Actual365Fixed(), jumps, jumpDates);
	euribor1m_curve->enableExtrapolation();
	

}

void get_spot_rates(ext::shared_ptr<YieldTermStructure> yieldCurve, DayCounter day_count, 
	Calendar calendar = UnitedStates(UnitedStates::FederalReserve), int months = 121) {
	vector<Real> spots = {};
	vector<Real> tenors = {};
	Date ref_date = yieldCurve->referenceDate();
	Date calc_date(ref_date);
	for (int month = 0; month< months; month++) {
		Real yrs = month / 12.0;
		Date d = calendar.advance(ref_date, Period(month, Months));
		Compounding compounding = Compounded;
		Frequency freq = Semiannual;
		InterestRate zeroRate = yieldCurve->zeroRate(yrs, compounding, freq);
		tenors.push_back(yrs);
		Real eq_rate = zeroRate.equivalentRate(day_count,compounding, freq, calc_date, d).rate();
		spots.push_back(100*eq_rate);
	}
	cout << "Maturities   Curve" << endl;
	for (int i = 116; i < tenors.size(); i++) {
		cout << tenors[i] << " " << spots[i] << endl;
	}



}
void yield_curve() {
	//the construction of treasury yield curve
	vector<Period> depo_maturities = { Period(6,Months), Period(12, Months) };
	vector<Rate> depo_rates = { 5.25, 5.5 };

	//bond rates
	vector<Rate> bond_rates = { 5.75, 6.0, 6.25, 6.5, 6.75, 6.80, 7.00, 7.1, 7.15, 7.2, 7.3, 7.35, 7.4, 7.5, 7.6, 7.6, 7.7, 7.8 };
	vector<Period> bond_maturities = {};
	for (int i = 3; i < 21; i++) {
		bond_maturities.push_back(Period(6 * i, Months));
	}

	vector<Period> maturities = depo_maturities;
	maturities.insert(maturities.end(), bond_maturities.begin(), bond_maturities.end());
	vector<Rate> rates = depo_rates;
	rates.insert(rates.end(), bond_rates.begin(), bond_rates.end());

	Date calc_date(15, January, 2015);
	Settings::instance().evaluationDate() = calc_date;

	Calendar calendar = UnitedStates(UnitedStates::SOFR);
	BusinessDayConvention business_convention = Unadjusted;
	DayCounter day_count = Thirty360(Thirty360::BondBasis);
	bool end_of_month = true;
	Natural settlement_days = 0;
	Real face_amount = 100;
	Period coupon_frequency(Semiannual);
	
	vector<ext::shared_ptr<RateHelper>> depo_helpers; 
	for (int i = 0; i < depo_rates.size(); i++) { 
		depo_helpers.push_back(ext::make_shared<DepositRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(depo_rates[i]/100.0)), 
			depo_maturities[i],  settlement_days, calendar, business_convention, end_of_month, day_count)); 
	}

	vector<ext::shared_ptr<RateHelper>> bond_helpers;
	for (int i = 0; i < bond_rates.size(); i++) {
		Date termination_date = calc_date + bond_maturities[i];   
		Schedule schedule(calc_date, termination_date, coupon_frequency, calendar, business_convention, business_convention, DateGeneration::Backward, end_of_month);
		auto bond_helper = ext::make_shared<FixedRateBondHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(face_amount)), settlement_days, face_amount,
			schedule, vector<Rate>(1,bond_rates[i] / 100.0), day_count, business_convention);
	//	auto test = bond_helper->quote()->value();
		bond_helpers.push_back(bond_helper);
		
	}

	vector<ext::shared_ptr<RateHelper>> rate_helpers(depo_helpers);
	rate_helpers.insert(rate_helpers.end(), bond_helpers.begin(), bond_helpers.end());
	
	auto yc_logcubicdiscount = ext::make_shared<PiecewiseYieldCurve<Discount, MonotonicLogCubic>>(calc_date, rate_helpers, day_count);
	yc_logcubicdiscount->enableExtrapolation();
	get_spot_rates(yc_logcubicdiscount, day_count);

	cout << "-----PiecewiseLinearZero-----" << endl;
	auto yc_linearzero = ext::make_shared<PiecewiseYieldCurve<ZeroYield, Linear>>(calc_date, rate_helpers, day_count);
	yc_linearzero->enableExtrapolation();
	get_spot_rates(yc_linearzero, day_count);

	cout << "-----PiecewiseCubicZero-----" << endl;
	auto yc_cubiczero = ext::make_shared<PiecewiseYieldCurve<ZeroYield, Cubic>>(calc_date, rate_helpers, day_count);
	yc_cubiczero->enableExtrapolation();
	get_spot_rates(yc_logcubicdiscount, day_count);

}

void implied_term_structures() {
	
	Date today(9, March, 2016);
	Settings::instance().evaluationDate() = today;

	vector<ext::shared_ptr<RateHelper>> helpers;

	vector<Period> tenor = { Period(6, Months), Period(2, Years), Period(5, Years), Period(10, Years), Period(15, Years) };
	vector<Rate> rates = { 0.201, 0.258, 0.464, 1.151, 1.588 };

	for (int i = 0; i < tenor.size(); i++) {
		auto helper = ext::make_shared<SwapRateHelper>(Handle<Quote>(ext::make_shared<SimpleQuote>(rates[i] / 100.0)), tenor[i], TARGET(), Annual,
			Unadjusted, Thirty360(Thirty360::BondBasis), ext::make_shared<Euribor6M>());
		helpers.push_back(helper);
	}
	auto curve = ext::make_shared<PiecewiseYieldCurve<ZeroYield, Linear>>(0, TARGET(), helpers, Actual360());

	//I'm using linear interpolation on the zero rates, which isn't great for actual use, However, the resulting jumps in the forward
	// rates will be useful as visual points of reference;note, for instance the jump around March 2018

	//The curve also implies an interest-rate curve in 1 year, 
	Date future_reference = today + Period(1, Years);
	auto implied_curve = ext::make_shared<ImpliedTermStructure>(Handle<YieldTermStructure>(curve),future_reference);

	vector<Date> dates;
	vector<Rate> rates_1, rates_2, rates_3;
	for (int i = 0; i < 6; i++) {
		Date date(future_reference + Period(i, Years));
		dates.push_back(date);
		rates_1.push_back(curve->forwardRate(future_reference, date,Actual360(),Continuous).rate());
		rates_2.push_back(implied_curve->zeroRate(date, Actual360(), Continuous).rate());
		cout << date << "   " << rates_1[i] << "   " << rates_2[i] << endl;
	}

	Settings::instance().evaluationDate() = future_reference; 
	cout << "------------date change---------------------" << endl;
	for (int i = 0; i < dates.size(); i++) {
		rates_3.push_back(implied_curve->zeroRate(dates[i],Actual360(), Continuous).rate());
		cout << dates[i] << "  " << rates_2[i] << "    " << rates_3[i] << endl;
	}
	//the reference date of the original curve was specified relative to the evaluation date, and when we moved it the curve moved, too.

	//after moving the evaluation date, the original and implied curve are exactly the same, and the spot rates return by the implied curve
	// are no longer forward rates, but the spot rates returned by the original curve

	cout << "----------moving the evaluation dates--------" << endl;
	rates_1.clear();
	rates_2.clear();
	for (int i = 0; i < dates.size(); i++) {
		rates_1.push_back(curve->forwardRate(future_reference, dates[i], Actual360(), Continuous).rate());
		rates_2.push_back(implied_curve->zeroRate(dates[i], Actual360(), Continuous).rate());
		cout << dates[i]<<"     "<<rates_1[i] << "   " << rates_2[i] << endl;
	}

	Settings::instance().evaluationDate() = today;
	

	vector<Date> node_dates;
	vector<Real> node_rates;
	auto nodes = curve->nodes();
	cout << endl;
	for (auto &node : nodes) {
		node_dates.push_back(node.first);
		node_rates.push_back(node.second);
		cout << node.first << " " << node.second << endl;
	}

	auto frozen_curve = ext::make_shared<InterpolatedZeroCurve<Linear>>(node_dates, node_rates, curve->dayCounter());
	implied_curve = ext::make_shared<ImpliedTermStructure>(Handle<YieldTermStructure>(frozen_curve), future_reference);
	Settings::instance().evaluationDate() = future_reference;

	cout << endl;
	for (int i = 0; i < dates.size(); i++) {
		Rate rate1 = frozen_curve->zeroRate(dates[i], Actual360(), Continuous).rate();
		Rate rate2 = frozen_curve->forwardRate(future_reference, dates[i], Actual360(), Continuous).rate();
		Rate rate3 = implied_curve->zeroRate(dates[i], Actual360(), Continuous).rate();
		cout << dates[i] << " " << rate1 << " " << rate2 << " " << rate3 << endl;
	}
}

void interest_rate_sensitivities() {
	Date today(8, March, 2016);
	Settings::instance().evaluationDate() = today;
	vector<ext::shared_ptr<SimpleQuote>> quotes;
	vector<ext::shared_ptr<RateHelper>> helpers;
	quotes.push_back(ext::make_shared<SimpleQuote>(0.312/100));
	helpers.push_back(ext::make_shared<DepositRateHelper>(Handle<Quote>(quotes[0]), Period(6, Months), 3, TARGET(), Following, false, Actual360()));

	vector<Real> rates = { 0.293, 0.272, 0.260, 0.256, 0.252, 0.248, 0.254, 0.261, 0.267, 0.279, 0.291, 0.303, 0.318, 0.335, 0.352, 0.371, 0.389, 0.409 };
	vector<int> months_to_start = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };

	for (int i = 0; i < rates.size(); i++) {
		quotes.push_back(ext::make_shared<SimpleQuote>(rates[i]/100));
		helpers.push_back(ext::make_shared<FraRateHelper>(Handle<Quote>(quotes[i+1]),months_to_start[i], ext::make_shared<Euribor6M>()));
	}

	rates = { 0.424, 0.576, 0.762, 0.954, 1.135, 1.303, 1.452, 1.584, 1.809, 2.037, 2.187, 2.234, 2.256, 2.295, 2.348, 2.421, 2.463 };

	// 期限数组
	vector<int> tenor = { 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30, 35, 40, 50, 60 };
	for (int i = 0; i < rates.size(); i++) {
		auto quote = ext::make_shared<SimpleQuote>(rates[i] / 100);
		quotes.push_back(quote);
		helpers.push_back(ext::make_shared<SwapRateHelper>(Handle<Quote>(quote), 
			Period(tenor[i], Years), TARGET(), Annual, Unadjusted, Thirty360(Thirty360::BondBasis), ext::make_shared<Euribor6M>()));
	}

	auto rate_curve = ext::make_shared<PiecewiseYieldCurve<Discount, MonotonicLogCubic>>(2, TARGET(), helpers, Actual365Fixed());
	RelinkableHandle<YieldTermStructure> curve_handle;
	curve_handle.linkTo(rate_curve);

	Schedule fixed_schedule(Date(8, April, 2016), Date(8, April, 2028), Period(1, Years), TARGET(), Following, Following, DateGeneration::Forward,false);
	Schedule floating_schedule(Date(8, April, 2016), Date(8, April, 2028), Period(1, Years), TARGET(), Following, Following, DateGeneration::Forward,false);
	auto index = ext::make_shared<Euribor6M>(curve_handle);

	VanillaSwap swap(VanillaSwap::Payer, 10000.0, fixed_schedule, 0.02, Thirty360(Thirty360::BondBasis), floating_schedule, index, 0.0 , Actual360());
	swap.setPricingEngine(ext::make_shared<DiscountingSwapEngine>(curve_handle));
	Real p0 = swap.NPV();
	cout << p0 << endl;

	auto bp = 1.0e-4;
	Real ref = quotes[0]->value();
	quotes[0]->setValue(ref+1*bp);
	cout << swap.NPV() << endl;
	quotes[0]->setValue(ref);

	for (auto& quote : quotes) {
		quote->setValue(quote->value() + 1 * bp);
	}
	cout << swap.NPV() << endl;
	for (auto& quote : quotes) {
		quote->setValue(quote->value() - 1 * bp);
	}
	// to shift all the zero rates upwards, you can keep the original curve as it is and add one basis point to all its zero rates by 
	// means of the ZeroSpreadedTermStructure class, To use it for pricing our swap, we'll store the original curve in a separate handle,
	// add the spread, and link the new curve to the handle we're using for forecasting.
	Handle<YieldTermStructure> base_curve(rate_curve);
	auto spread = ext::make_shared<SimpleQuote>(1 * bp);
	curve_handle.linkTo(ext::make_shared<ZeroSpreadedTermStructure>(base_curve, Handle<Quote>(spread)));
	cout << swap.NPV() << endl;

	Date spot = rate_curve->referenceDate();
	vector<Date> dates;
	vector<Handle<Quote>> spreads;
	for (int i = 0; i < 21; i++) {
		dates.push_back(spot+Period(i, Years));
		spreads.push_back(Handle<Quote>(ext::make_shared<SimpleQuote>((i-7)*bp)));
	}
	curve_handle.linkTo(ext::make_shared<InterpolatedPiecewiseZeroSpreadedTermStructure<Linear>>(base_curve,spreads,dates));
	cout << swap.NPV() << endl;

	curve_handle.linkTo(rate_curve);
	cout << swap.NPV() << endl;

}

void forward_rate_curves() {
	Date today(24, August, 2015);
	Settings::instance().evaluationDate() = today;

	vector<int> tenor = { 1, 2, 3, 5, 10, 20 };
	vector<Rate> rates = {0.01, 0.03, 0.02, 0.025, 0.035, 0.05, 0.04};
	vector<Date> dates = { today };
	for (int i = 0; i < tenor.size(); i++) 
		dates.push_back(today + Period(tenor[i], Years));
	ForwardCurve curve(dates, rates, Actual360());
	auto nodes = curve.nodes();
	Date d = today + Period(4, Years);
	cout << endl << d << " " << curve.forwardRate(d, d, curve.dayCounter(), Continuous) << endl << endl;

	dates.clear(); 
	rates.clear(); 
	vector<Real> expected; 
	for (int i = 0; i < nodes.size(); i++) { 
		dates.push_back(nodes[i].first); 
		rates.push_back(curve.forwardRate(dates[i], dates[i], curve.dayCounter(), Continuous).rate());
		expected.push_back(nodes[i].second);
		cout << dates[i] << "  " << expected[i] << "  " << rates[i] << endl;
	}

}