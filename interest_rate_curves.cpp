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
		synth_helpers.push_back(ext::make_shared<SwapRateHelper>(Handle<Quote>()));
	}




	






}

void euribor_curve2() {
	

}