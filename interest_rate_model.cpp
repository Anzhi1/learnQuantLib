#include "interest_rate_model.h"
#include <ql/quantlib.hpp>
#include <iostream>
#include <cmath>

using namespace std;
using namespace QuantLib;

typedef RandomSequenceGenerator<CLGaussianRng<MersenneTwisterUniformRng>> GSG;
typedef RandomSequenceGenerator<InverseCumulativeRng<MersenneTwisterUniformRng, InverseCumulativeNormal>> ISG;
typedef RandomSequenceGenerator<InverseCumulativeRng<SobolRsg, InverseCumulativeNormal>> SSG;

void hull_white_model() {
	Real sigma = 0.1;
	Real a = 0.1;
	Natural timestamp = 360;
	Natural length = 30;//in years
	Rate forward_rate = 0.05;
	DayCounter day_count = Thirty360(Thirty360::BondBasis);
	Date todays_date(15, January, 2015);

	Settings::instance().evaluationDate() = todays_date;
	auto spot_curve = ext::make_shared<FlatForward>(todays_date, Handle<Quote>(ext::make_shared<SimpleQuote>(forward_rate)), day_count);
	Handle<YieldTermStructure> spot_curve_handle(spot_curve);
	ext::shared_ptr<StochasticProcess> hw_process = ext::make_shared<HullWhiteProcess>(spot_curve_handle, a, sigma);
	//create mersenne twister uniform  random generator
	unsigned long seed = 28749;
	MersenneTwisterUniformRng generator(seed);
	//GaussianRandomSequenceGenerator

	// template <class URNG, class IC>
	// struct GenericPseudoRandom
	// typedef URNG urng_type
	// typedef InverseCumulativeRng<urng_type,IC> rng_type;
	// typedef RandomSequenceGenerator<urng_type> ursg_type;
	// typedef InverseCumulativeRsg<ursg_type, IC> rsg_type;

	// GaussianPathGenerator-->PathGenerator<GaussianRandomSequenceGenerator>-->
	// PathGenerator<InverseCumulativeRsg<RandomSequenceGenerator<UniformRandomGenerator>>>
	// 
	// GaussianRandomSequenceGenerator-->rsg_type-->InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>,InverseCumulativeNormal>
	// UniformRandomSequenceGenerator-->ursg_type
	// UniformRandomGenerator--> urng_type
	//create gaussian generator by using central limit transformation method
	//CLGaussianRng<MersenneTwisterUniformRng> gaussianGenerator(generator);
	InverseCumulativeRng<MersenneTwisterUniformRng, InverseCumulativeNormal> mygenerator(generator);

	//define maturity, number of steps per path and create gaussian sequence generator
	//GSG gaussianSequenceGenerator(timestamp, gaussianGenerator);
	ISG mySequenceGenerator(timestamp, mygenerator);

	//PathGenerator<GSG> pathGenerator(hw_process, length, timestamp, gaussianSequenceGenerator, false);
	PathGenerator<ISG> path_generator(hw_process, length, timestamp, mySequenceGenerator, false);

	//create matrix 
	Size nColumns = 10;
	Matrix paths(timestamp + 1, nColumns);
	for (unsigned int i = 0; i != paths.columns(); i++) {
		//request a new stochastic path from path generator
		Sample<Path> path = path_generator.next();
		//		cout << path.value.length() << endl;
		//		cout << path.weight << endl;

				//save generated path into container
		for (unsigned int j = 0; j != path.value.length(); j++) {
			paths[j][i] = path.value.at(j);
			//		cout << paths[j][i] << " ";
		}
		//	cout << endl;
	}
	cout << paths.rows() << " " << paths.columns() << endl;
	for (int i = 0; i < paths.rows(); i++) {
		for (int j = 0; j < paths.columns(); j++) {
			cout << to_string(paths[i][j]) << " ";
		}
		cout << endl;
	}
}

// returns a path generator
template<typename SSG, typename ISG>
PathGenerator<SSG> get_path_generator(Natural timestep, ext::shared_ptr<StochasticProcess>hw_process, Natural length, bool low_discrepancy, bool brownian_bridge = true) {
	if (low_discrepancy) {
		SobolRsg sobolGen(1);
		InverseCumulativeRng<SobolRsg, InverseCumulativeNormal> generator(sobolGen);
		SSG rng(timestep, generator);
		PathGenerator<SSG> seq(hw_process, length, timestep, rng, brownian_bridge);
	}
	else {
		MersenneTwisterUniformRng ugenerator(1);
		InverseCumulativeRng<MersenneTwisterUniformRng, InverseCumulativeNormal> generator(ugenerator);
		ISG rng(timestep, generator);
		PathGenerator<ISG> seq(hw_process, length, timestep, rng, brownian_bridge);
	}


}
void Monte_Carlo_Hull_White() {

	Date todays_date(15, January, 2015);
	Settings::instance().evaluationDate() = todays_date;

}

typedef std::tuple<double, double, double> CalibrationData;

vector<ext::shared_ptr<CalibrationHelper>> create_swaption_helpers(vector<CalibrationData> data,
	ext::shared_ptr<Euribor1Y> index,
	Handle<YieldTermStructure> term_structure,
	ext::shared_ptr<PricingEngine> engine) {

	vector<ext::shared_ptr<CalibrationHelper>> swaptions;
	Period fixed_leg_tenor(1, Years);
	DayCounter fixed_leg_daycounter = Actual360();
	DayCounter floating_leg_daycounter = Actual360();

	for (const auto& tuple : data) {
		double start = get<0>(tuple);
		double length = get<1>(tuple);
		double volatility = get<2>(tuple);

		Handle<Quote> vol_handle(ext::make_shared<SimpleQuote>(volatility));
		auto helper = ext::make_shared<SwaptionHelper>(Period(start, Years), Period(length, Years),
			vol_handle, index, fixed_leg_tenor, fixed_leg_daycounter, floating_leg_daycounter, term_structure);
		helper->setPricingEngine(engine);
		swaptions.push_back(helper);
	}
	return swaptions;

}
vector<ext::shared_ptr<CalibrationHelper>> create_swaption_helpers_normal(vector<CalibrationData> data,
	ext::shared_ptr<Euribor1Y> index,
	Handle<YieldTermStructure> term_structure,
	ext::shared_ptr<PricingEngine> engine) {
	vector<ext::shared_ptr<CalibrationHelper>> swaptions;
	Period fixed_leg_tenor(1, Years);
	DayCounter fixed_leg_daycounter = Actual360();
	DayCounter floating_leg_daycounter = Actual360();

	for (const auto& tuple : data) {
		double start = get<0>(tuple);
		double length = get<1>(tuple);
		double volatility = get<2>(tuple);

		Handle<Quote> vol_handle(ext::make_shared<SimpleQuote>(volatility));
		auto helper = ext::make_shared<SwaptionHelper>(Period(start, Years), Period(length, Years),
			vol_handle, index, fixed_leg_tenor, fixed_leg_daycounter, floating_leg_daycounter, term_structure, BlackCalibrationHelper::RelativePriceError, Null<Real>(), 1.0, Normal);
		helper->setPricingEngine(engine);
		swaptions.push_back(helper);
	}

	return swaptions;
}
typedef std::tuple<double, double, double> CalibrationData;
void calibration_report(vector<ext::shared_ptr<SwaptionHelper>>& swaptions, vector<CalibrationData> data) {
	Real cum_err = 0.0, cum_err2 = 0.0;
	string s = "      ";
	cout << "Model Price" << s << "Market Price" << s << "Implied Vol" << s << "Market Vol" << "Rel Error Price" << s << "Rel Error Vols" << endl;
	for (int i = 0; i < swaptions.size(); i++) {
		auto a = swaptions[i];
		Real model_price = swaptions[i]->modelValue();
		Real market_vol = get<2>(data[i]);
		Real black_price = swaptions[i]->blackPrice(market_vol);
		Real rel_error = model_price / black_price - 1.0;
		Real implied_vol = swaptions[i]->impliedVolatility(model_price, 1e-5, 50, 0.0, 0.50);
		Real rel_error2 = implied_vol / market_vol - 1.0;
		cum_err += rel_error * rel_error;
		cum_err2 += rel_error2 * rel_error2;
		cout << model_price << s << black_price << s << implied_vol << s << market_vol << s << rel_error << s << rel_error2 << endl;
	}

	cout << "Cumulative Error Price: " << sqrt(cum_err) << endl;
	cout << "Cumulative Error Vols: " << sqrt(cum_err2) << endl;


}
//examples of calibrating the interest rate models to the swaption volatilities. we looked at setting up different interest rate models 
//and discussed both lognormal and normal volatilities
void short_interest_rate_model_calibration() {
	Date today(15, February, 2002);
	Date settlement(19, February, 2002);
	Settings::instance().evaluationDate() = today;
	Handle<YieldTermStructure> term_structure(ext::make_shared<FlatForward>(settlement, 0.04875825, Actual365Fixed()));
	auto index = ext::make_shared<Euribor1Y>(term_structure);

	//CalibrationData
	vector<CalibrationData> data = {
		make_tuple(1, 5, 0.1148),
		make_tuple(2, 4, 0.1108),
		make_tuple(3, 3, 0.1070),
		make_tuple(4, 2, 0.1021),
		make_tuple(5, 1, 0.1000)
	};

	//Hull-White 1 Factor Model
	auto model = ext::make_shared<HullWhite>(term_structure);
	auto engine = ext::make_shared<JamshidianSwaptionEngine>(model);
	vector<ext::shared_ptr<CalibrationHelper>> swaptions = create_swaption_helpers(data, index, term_structure, engine);
	LevenbergMarquardt optimization_method(1.0e-8, 1.0e-8, 1.0e-8);
	EndCriteria end_criteria(10000, 100, 1e-6, 1e-8, 1e-8);
	model->calibrate(swaptions, optimization_method, end_criteria);
	Array result = model->params();
	cout << result << endl;
	//  you cannot use the below method directly for some stupid reason about c++ , if u have better idea, tell me please. 
	//	calibration_report(swaptions,data);

		//Calibrating Volatility With Fixed Reversion¡ª¡ªperfrom calibration with constraints 
	auto constrained_model = ext::make_shared<HullWhite>(term_structure, 0.05, 0.001);
	engine = ext::make_shared<JamshidianSwaptionEngine>(constrained_model);
	swaptions = create_swaption_helpers(data, index, term_structure, engine);
	constrained_model->calibrate(swaptions, optimization_method, end_criteria, NoConstraint(), {}, vector<bool>{ true, false });
	result = constrained_model->params();
	cout << "constrained model calibrate result: " << result << endl;


	//Black Karasinski Model
	//Black Karasinski Model is not an affine model, and hence we cannot use the JamshidianSwaptionEngine. In order to calibrate, we use the
	//TreeSwaptionEngine which will work with all short rate models.
	auto bk_model = ext::make_shared<BlackKarasinski>(term_structure);
	auto bk_engine = ext::make_shared<TreeSwaptionEngine>(bk_model, 100);
	swaptions = create_swaption_helpers(data, index, term_structure, bk_engine);
	end_criteria = EndCriteria(10000, 100, 1e-6, 1e-8, 1e-8);
	bk_model->calibrate(swaptions, optimization_method, end_criteria);
	cout << "Black Karasinski Model: " << bk_model->params() << endl;


	//G2++ Model a calibration example of the 2-factor G2++ model
	auto g2_model = ext::make_shared<G2>(term_structure);
	auto g2_engine = ext::make_shared<TreeSwaptionEngine>(g2_model, 25);
	//	auto g2_engine = ext::make_shared<G2SwaptionEngine>(g2_model, 10, 400);
	//	auto g2_engine = ext::make_shared<TreeSwaptionEngine>(g2_model);
	end_criteria = EndCriteria(1000, 100, 1e-6, 1e-8, 1e-8);
	swaptions = create_swaption_helpers(data, index, term_structure, g2_engine);
	g2_model->calibrate(swaptions, optimization_method, end_criteria);
	cout << "G2 Model: " << g2_model->params() << endl;


	//Calibrating to Normal Volatilities
	swaptions = create_swaption_helpers_normal(data, index, term_structure, engine);
	end_criteria = EndCriteria(10000, 100, 1e-6, 1e-8, 1e-8);
	try {
		model->calibrate(swaptions, optimization_method, end_criteria);
		cout << "Normal Volatilities Hull-White model: " << model->params() << endl;
	}
	catch (exception& e) {
		cout << e.what() << endl;
	}

}


void par_versus_index_coupons() {

	Date today(7, January, 2013);
	Settings::instance().evaluationDate() = today;
	vector<Date> dates = {
		Date(7, January, 2013),
		Date(8, April, 2013),
		Date(8, July, 2013),
		Date(7, January, 2014),
		Date(7, July, 2014)
	};
	vector<double> forwards = {
		0.03613672438543303,
		0.03613672438543303,
		0.033849133719219514,
		0.03573931373272106,
		0.03445303757052511
	};
	auto libor_curve = ext::make_shared<ForwardCurve>(dates, forwards, Actual365Fixed());
	auto index = ext::make_shared<GBPLibor>(Period(6, Months), Handle<YieldTermStructure>(libor_curve));
	Calendar calendar = index->fixingCalendar();
	Real nominal = 1000000;
	Natural length = 1;
	Date maturity = calendar.advance(today, length, Years);
	BusinessDayConvention adjustment = index->businessDayConvention();
	Schedule schedule(today, maturity, index->tenor(), calendar, adjustment, adjustment, DateGeneration::Backward, false);
	IborLeg floating_leg(schedule, index);
	floating_leg.withPaymentDayCounter(index->dayCounter()).withNotionals(vector<Real>{nominal});

	dates = schedule.dates();

	vector<Date> fixing_date(dates.begin(), dates.end() - 1);
	vector<Real> index_fixing;
	for (int i = 0; i < fixing_date.size(); i++) {
		index_fixing.push_back(index->fixing(fixing_date[i]));
	}
	vector<Date> start_date(dates.begin(), dates.end() - 1);
	vector<Date> end_date(dates.begin() + 1, dates.end());
	vector<Real> days;
	vector<Real> accrual_period;
	vector<Real> amount;
	for (int i = 0; i < start_date.size(); i++) {
		auto diff = end_date[i] - start_date[i];
		days.push_back(diff);
		accrual_period.push_back(diff / 365.0);
		amount.push_back(index_fixing[i] * nominal * accrual_period[i]);
	}

	for (int i = 0; i < fixing_date.size(); i++) {
		cout << fixing_date[i] << "    " << index_fixing[i] << "    " << start_date[i] << "    " << end_date[i] << "    " << days[i] << "    " << accrual_period[i]
			<< "    " << amount[i] << endl;
	}
	auto cashflow = floating_leg.operator QuantLib::Leg();


	for (int i = 0; i < cashflow.size(); i++) {
		cout << cashflow[i]->amount() << "   " << ext::dynamic_pointer_cast<Coupon>(cashflow[i])->rate() << endl;
	}
	auto coupon = ext::dynamic_pointer_cast<FloatingRateCoupon>(cashflow[1]);
	cout << coupon->fixingDate() << "    " << index->fixing(coupon->fixingDate()) << endl;

	//the fixing is also consistent with what we can forecast from the LIBOR curve, given the start and end date of the underlying tenor:

	Date startDate = index->valueDate(coupon->fixingDate());
	Date endDate = index->maturityDate(startDate);

	cout << startDate << "   " << endDate << endl;

	cout << libor_curve->forwardRate(startDate, endDate, coupon->dayCounter(), Simple) << endl;

	//for historical reasons, the coupon is calculated at par; that is, the floating rate is calculated over the duration of the coupon, Due to the constraints of the schedule,
	//the end of the coupon doesn't correspond to the end of the LIBOR tenor
	Date couponStart(coupon->accrualStartDate());
	Date couponEnd(coupon->accrualEndDate());

	cout << couponStart << "   " << couponEnd << endl;

	cout << libor_curve->forwardRate(couponStart, couponEnd, coupon->dayCounter(), Simple) << endl;
	cout << coupon->rate() << endl;
	cout << coupon->rate() * nominal * coupon->accrualPeriod() << endl;
	cout << coupon->amount() << endl;
}

void caps_and_floors() {

	Date calc_date(14, June, 2016);
	Settings::instance().evaluationDate() = calc_date;
	vector<Date> dates = { Date(14,June,2016), Date(14,September,2016), Date(14,December,2016), Date(14,June,2017),  Date(14,June,2019),
		Date(14,June,2021), Date(15,June,2026), Date(16,June,2031), Date(16,June,2036), Date(14,June,2046) };
	vector<Real> yields = { 0.000000, 0.006616, 0.007049, 0.007795, 0.009599, 0.011203, 0.015068, 0.017583, 0.018998, 0.020080 };
	DayCounter day_count = ActualActual(ActualActual::ISDA);
	Calendar calendar = UnitedStates(UnitedStates::GovernmentBond);
	auto interpolation = Linear();
	Compounding compounding = Compounded; 
	Frequency compounding_frequency = Annual;
	auto term_structure = ext::make_shared<InterpolatedZeroCurve<Linear>>(dates, yields, day_count, calendar, interpolation, compounding, compounding_frequency);
	auto ts_handle = Handle<YieldTermStructure>(term_structure);

	//as a next step,lets construct the cap itself, In order to do that, we start by constructing the Schedule object to project the cash flows
	Date start_date(14, June, 2016);
	Date end_date(14, June, 2026);
	Period period(3, Months);
	BusinessDayConvention buss_convention = ModifiedFollowing;
	DateGeneration::Rule rule = DateGeneration::Forward;	
	bool end_of_month = false;

	Schedule schedule(start_date, end_date, period, calendar, buss_convention, buss_convention, rule, end_of_month);
	auto ibor_index = ext::make_shared<USDLibor>(Period(3, Months), ts_handle);
	ibor_index->addFixing(Date(10,June,2016), 0.0065560);
	IborLeg ibor_leg(schedule, ibor_index);
	ibor_leg.withNotionals(vector<Real>{1000000});
	
	Real strike = 0.02;
	Cap cap(ibor_leg, vector<Real>{strike});
	Handle<Quote> vol(ext::make_shared<SimpleQuote>(0.547295));
	auto engine = ext::make_shared<BlackCapFloorEngine>(ts_handle, vol);
	cap.setPricingEngine(engine);
	cout<<setw(12) << setprecision(12) << cap.NPV() << endl;

	
	//Using Volatility Surfaces
	vector<Real> strikes = { 0.01, 0.015, 0.02 };
	vector<int> temp = { 1,2,3,4,5,6,7,8,9,10,12 };
	vector<Period> expires;
	for (auto& t : temp) {
		expires.push_back(Period(t,Years));
	}
	Matrix vols(expires.size(), strikes.size());
	vector<vector<Real>> data = {
		{ 47.27, 55.47, 64.07, 70.14, 72.13, 69.41, 72.15, 67.28, 66.08, 68.64, 65.83 },
		{ 46.65, 54.15, 61.47, 65.53, 66.28, 62.83, 64.42, 60.05, 58.71, 60.35, 55.91 },
		{ 46.6,  52.65, 59.32, 62.05, 62.0,  58.09, 59.03, 55.0,  53.59, 54.74, 49.54 }
	};
	for (int i = 0; i < vols.rows(); i++) {
		for (int j = 0; j < vols.columns(); j++) {
			vols[i][j] = data[j][i] / 100.0;
		}
	}

	BusinessDayConvention bdc = ModifiedFollowing;
	day_count = Actual365Fixed();
	Natural settlement_days = 2;
	auto capfloor_vol = ext::make_shared<CapFloorTermVolSurface>(settlement_days, calendar, bdc, expires, strikes, vols, day_count); 
	
	OptionletStripper1 optionlet_surf(capfloor_vol, ibor_index, {}, 1e-6, 100, ts_handle);
	Handle<OptionletVolatilityStructure> ovs_handle(ext::make_shared<StrippedOptionletAdapter>(ext::make_shared<OptionletStripper1>(optionlet_surf)));

	auto engine2 = ext::make_shared<BlackCapFloorEngine>(ts_handle,ovs_handle);
	cap.setPricingEngine(engine2);
	cout << cap.NPV() << endl;
	cout << "infer the implied volatility for the cap at its NPV: " << cap.impliedVolatility(cap.NPV(), ts_handle, 0.4) << endl;






}