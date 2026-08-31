// Unit tests for TIM::Time, TIM::Date, and TIM::Calendar

#include <chrono>

#include <gtest/gtest.h>

#include "core/tim_time.hpp"

namespace {

using namespace std::chrono_literals;
using TIM::Calendar;
using TIM::Date;
using TIM::Duration;
using TIM::Time;

constexpr Calendar kCalendars[] = {Calendar::ThirtyDayMonths, Calendar::Julian,
                                   Calendar::Gregorian, Calendar::NoLeap};

// The base date of every calendar is the start of year 1, i.e. a zero span.
TEST(Time, BaseDateIsZero) {
    for (const Calendar cal : kCalendars) {
        EXPECT_EQ(Date{}.to_time(cal), Time{});
        EXPECT_EQ(Time{}.to_date(cal), Date{});
    }
}

// Date -> Time -> Date is the identity.
TEST(Time, DateRoundTrips) {
    const Date dates[] = {{1, 1, 1}, {1, 12, 28, 23, 59, 59}, {1900, 3, 1},
                          {2000, 2, 28, 6, 30, 15}, {2024, 12, 28}, {4000, 7, 4}};
    for (const Calendar cal : kCalendars)
        for (const Date& d : dates)
            EXPECT_EQ(d.to_time(cal).to_date(cal), d) << d.to_string();
}

// Round trips at the ends of a month, where the calendars disagree on which
// dates exist at all.
TEST(Time, EndOfMonthRoundTrips) {
    const Date gregorian[] = {{2024, 2, 29}, {2023, 2, 28}, {2024, 12, 31}, {1900, 2, 28}};
    for (const Date& d : gregorian)
        EXPECT_EQ(d.to_time(Calendar::Gregorian).to_date(Calendar::Gregorian), d) << d.to_string();
    const Date thirty{2024, 2, 30};
    EXPECT_EQ(thirty.to_time(Calendar::ThirtyDayMonths).to_date(Calendar::ThirtyDayMonths), thirty);
    const Date julian{1900, 2, 29};  // leap in Julian, not in Gregorian
    EXPECT_EQ(julian.to_time(Calendar::Julian).to_date(Calendar::Julian), julian);
}

// Each calendar's leap rule, read through the month lengths it implies.
TEST(Calendar, MonthLengthsFollowTheLeapRule) {
    EXPECT_EQ(TIM::days_in_month(Calendar::Gregorian, 2024, 2), 29);
    EXPECT_EQ(TIM::days_in_month(Calendar::Gregorian, 2023, 2), 28);
    EXPECT_EQ(TIM::days_in_month(Calendar::Gregorian, 1900, 2), 28);  // 100 not 400
    EXPECT_EQ(TIM::days_in_month(Calendar::Gregorian, 2000, 2), 29);  // 400
    EXPECT_EQ(TIM::days_in_month(Calendar::Julian, 1900, 2), 29);     // every 4th
    EXPECT_EQ(TIM::days_in_month(Calendar::NoLeap, 2024, 2), 28);
    EXPECT_EQ(TIM::days_in_month(Calendar::ThirtyDayMonths, 2024, 2), 30);
    EXPECT_EQ(TIM::days_in_month(Calendar::Gregorian, 2024, 12), 31);
}

// Julian and Gregorian have drifted apart by year 1900, so the same date is a
// different point on the axis in each.
TEST(Calendar, JulianAndGregorianDiverge) {
    const Date d{1900, 3, 1};
    EXPECT_NE(d.to_time(Calendar::Julian), d.to_time(Calendar::Gregorian));
    const Date base{1, 1, 1};
    EXPECT_EQ(base.to_time(Calendar::Julian), base.to_time(Calendar::Gregorian));
}

// Exact-unit arithmetic is Duration arithmetic: it takes no calendar, and the
// compiler converts coarser units for free.
TEST(Time, DurationArithmeticNeedsNoCalendar) {
    const Time t = Date{2024, 2, 28, 12, 0, 0}.to_time(Calendar::Gregorian);
    EXPECT_EQ((t + 24h).to_date(Calendar::Gregorian).day, 29);  // into the leap day
    EXPECT_EQ((t + std::chrono::days{2}).to_date(Calendar::Gregorian).month, 3);
    EXPECT_EQ(t + 1h - 1h, t);
    EXPECT_LT(t, t + 1s);
}

// Differencing two points yields an exact span, the building block for a
// "seconds since <base>" time axis.
TEST(Time, DifferenceIsExact) {
    const Time a = Date{2024, 1, 1}.to_time(Calendar::Gregorian);
    const Time b = Date{2024, 1, 2}.to_time(Calendar::Gregorian);
    EXPECT_EQ(b - a, 86400s);
    EXPECT_EQ(a - b, -86400s);
    EXPECT_EQ(b - b, Duration::zero());
}

// Calendar increments keep the day of month, and borrow across a year end.
TEST(Time, CalendarIncrements) {
    const Calendar cal = Calendar::Gregorian;
    const Time t = Date{2024, 1, 15}.to_time(cal);
    EXPECT_EQ(t.add_months(cal, 1).to_date(cal), (Date{2024, 2, 15}));
    EXPECT_EQ(t.add_months(cal, 12).to_date(cal), (Date{2025, 1, 15}));
    EXPECT_EQ(t.add_months(cal, -1).to_date(cal), (Date{2023, 12, 15}));
    EXPECT_EQ(t.add_years(cal, -1).to_date(cal), (Date{2023, 1, 15}));
    // A leap day survives a jump of four years, which is why add_years keeps
    // the day rather than clamping it.
    const Time leap = Date{2024, 2, 29}.to_time(cal);
    EXPECT_EQ(leap.add_years(cal, 4).to_date(cal), (Date{2028, 2, 29}));
}

// The enumerator values are the FMS time_manager parameters, so the bridge can
// pass get_calendar_type() straight through.
TEST(Calendar, MatchesFmsParameterValues) {
    EXPECT_EQ(TIM::calendar_from_fms(0), Calendar::NoCalendar);
    EXPECT_EQ(TIM::calendar_from_fms(1), Calendar::ThirtyDayMonths);
    EXPECT_EQ(TIM::calendar_from_fms(2), Calendar::Julian);
    EXPECT_EQ(TIM::calendar_from_fms(3), Calendar::Gregorian);
    EXPECT_EQ(TIM::calendar_from_fms(4), Calendar::NoLeap);
}

}  // namespace
