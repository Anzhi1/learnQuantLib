#include <iostream>
#include <ql/quantlib.hpp>
#include "basics.h"

using namespace std;
using namespace QuantLib;

void date_class() {

    //1.Date Class
    Date date(31, March, 2015);
    cout << date << endl;

    cout << date.month() << "-" << date.dayOfMonth() << "-" << date.year() << endl;

    cout << date.weekday() << endl;

    cout << (date.weekday() == Tuesday);

    cout << "add a day     :" << date + 1 << endl;
    cout << "subtract a day:" << date - 1 << endl;
    cout << "add a week    :" << date + Period(1, Weeks) << endl;
    cout << "add a month   :" << date - Period(1, Months) << endl;
    cout << "add a year    :" << date + Period(1, Years) << endl;

    cout << (date == Date(31, March, 2015)) << endl;
    cout << (date > Date(31, March, 2015)) << endl;
    cout << (date < Date(1, April, 2015)) << endl;
    cout << (date != Date(1, April, 2015));
}

void calendar_class() {
    Date date(31, March, 2015);
    // to do
    Calendar us_calendar = UnitedStates(UnitedStates::Settlement);
    Calendar italy_calendar = Italy(Italy::Settlement);

    Period period(60, Days);
    Date raw_date = date + period;
    Date us_date = us_calendar.advance(date, period);
    Date italy_date = italy_calendar.advance(date, period);

    cout << "add 60 days:" << raw_date << endl;
    cout << "add 60 business days in US:" << us_date << endl;
    cout << "add 60 business days in Italy:" << italy_date << endl;


    auto us_busdays = us_calendar.businessDaysBetween(date, us_date);
    auto italy_busdays = italy_calendar.businessDaysBetween(date, italy_date);

    cout << "Buisness days US: " << us_busdays << endl;
    cout << "Business days Italy:" << italy_busdays << endl;

    //joint calendar 

    JointCalendar joint_calendar(us_calendar, italy_calendar);

    Date joint_date = joint_calendar.advance(date, period);
    auto joint_busdays = joint_calendar.businessDaysBetween(date, joint_date);

    cout << "add 60 business days in US-Italy:" << joint_date << endl;
    cout << "Business days US-Italy:" << joint_busdays << endl;

}

void print_schedule(Schedule schedule) {
    vector<Date> dates = schedule.dates();
    for (int i = 0; i < dates.size(); i++) {
        cout << i + 1 << " " << dates[i] << endl;
    }
    
}

void schedule_class() {
    Date effective_date(1, January, 2015);
    Date termination_date(1, January, 2016);
    Period tenor(Monthly);
    Calendar calendar = UnitedStates(UnitedStates::Settlement);
    BusinessDayConvention business_convention = Following;
    BusinessDayConvention termination_business_convention = Following;
    bool end_of_month = false;

    Schedule schedule(effective_date, termination_date, tenor, calendar, business_convention, termination_business_convention, DateGeneration::Backward, end_of_month);
    print_schedule(schedule);
    cout << "----------short stub in the front----------" << endl;
    //short stub in the front
    Date first_date(15, January, 2015);
    Schedule schedule2(effective_date, termination_date, tenor, calendar, business_convention, termination_business_convention, DateGeneration::Backward, end_of_month, first_date);
    print_schedule(schedule2);

    //short stub at the back
    cout << "----------short stub at the back----------" << endl;
    Date penultimate_date(15, December, 2015);
    Schedule schedule3(effective_date, termination_date, tenor, calendar, business_convention, termination_business_convention, DateGeneration::Forward, end_of_month, Date(), penultimate_date);
    print_schedule(schedule3);

    //long stub in the front
    cout << "----------long stub in the front----------" << endl;
    first_date = Date(1, February, 2015);
    effective_date = Date(15, December, 2014);
    termination_date = Date(1, January, 2016);
    Schedule schedule4(effective_date, termination_date, tenor, calendar, business_convention, termination_business_convention, DateGeneration::Backward, end_of_month, first_date);
    print_schedule(schedule4);

    //long stub at the back
    cout << "----------long stub at the back----------" << endl;
    effective_date = Date(15, January, 2015);
    termination_date = Date(15, January, 2016);
    penultimate_date = Date(1,December,2015);
    Schedule schedule5(effective_date, termination_date, tenor, calendar, business_convention, termination_business_convention, DateGeneration::Forward, end_of_month, Date(),penultimate_date);
    print_schedule(schedule5);


    //build schedule by a list of Date
    cout << "----------build schedule by a list of Date----------" << endl;
    vector<Date> dates;
    dates.push_back(Date(2, January, 2015));    dates.push_back(Date(2, February, 2015));
    dates.push_back(Date(2, March, 2015));      dates.push_back(Date(1, April, 2015));
    dates.push_back(Date(1, May, 2015));        dates.push_back(Date(1, June, 2015));
    dates.push_back(Date(1, July, 2015));       dates.push_back(Date(3, August, 2015));
    dates.push_back(Date(1, September, 2015));  dates.push_back(Date(1, October, 2015));
    dates.push_back(Date(2, November, 2015));   dates.push_back(Date(2, December, 2015));
    dates.push_back(Date(4, January, 2016));

    Schedule schedule6(dates);
    print_schedule(schedule6);

}

void interest_rate() {
    Real annual_rate = 0.05;
    DayCounter dc = ActualActual(ActualActual::ISDA);
    Compounding compound_type = Compounded;
    Frequency frequency = Annual;
    InterestRate interest_rate(annual_rate, dc, compound_type, frequency);

    Real t = 2.0;
    //compoundFactor，复利因子，若为simple,则为1+rt
    cout << interest_rate.compoundFactor(t)<<endl;
    cout << (1 + annual_rate) * (1.0 + annual_rate)<<endl;

    cout << setprecision(10);
    cout<< interest_rate.discountFactor(t) << endl;
    cout << 1.0 / interest_rate.compoundFactor(t) << endl;

    Frequency new_frequency = Semiannual;
    InterestRate new_interest_rate = interest_rate.equivalentRate(compound_type, new_frequency, t);
    Real new_annual_rate = new_interest_rate.rate();
    cout << new_annual_rate << endl;

    cout << interest_rate.discountFactor(t) << endl;
    cout << new_interest_rate.discountFactor(t) << endl;

}

void instruments() {
    Date today(7, March, 2014);
    Settings::instance().evaluationDate() = today;



}


    



    
    