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
	auto spot_curve = ext::make_shared<FlatForward>(todays_date, Handle<Quote>(ext::make_shared<SimpleQuote>(forward_rate)),day_count);
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
			cout << to_string(paths[i][j])<< " ";
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
		SSG rng(timestep,generator);
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
		auto helper = ext::make_shared<SwaptionHelper>(Period(start,Years), Period(length, Years),
			vol_handle, index, fixed_leg_tenor, fixed_leg_daycounter, floating_leg_daycounter, term_structure);
		helper->setPricingEngine(engine);
		swaptions.push_back(helper);
	}
	return swaptions;

}
typedef std::tuple<double, double, double> CalibrationData;
void calibration_report(vector<ext::shared_ptr<SwaptionHelper>> &swaptions, vector<CalibrationData> data) {
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

	auto model = ext::make_shared<HullWhite>(term_structure);
	auto engine = ext::make_shared<JamshidianSwaptionEngine>(model);
	vector<ext::shared_ptr<CalibrationHelper>> swaptions = create_swaption_helpers(data, index, term_structure, engine);
	LevenbergMarquardt optimization_method(1.0e-8, 1.0e-8, 1.0e-8);
	EndCriteria end_criteria(10000, 100, 1e-6, 1e-8, 1e-8);
	model->calibrate(swaptions, optimization_method, end_criteria);
	auto result = model->params();
	cout << result << endl;
//  you cannot use the below method directly for some stupid reason about c++ , if u have better idea, tell me please. 
//	calibration_report(swaptions,data);

	//Calibrating Volatility With Fixed Reversion¡ª¡ªperfrom calibration with constraints 
	auto constrained_model = ext::make_shared<HullWhite>(term_structure, 0.05, 0.001);





}