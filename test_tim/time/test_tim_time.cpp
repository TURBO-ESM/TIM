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

// add_months keeps the day of month, so an increment landing on a day the
// target month does not have is fatal.
TEST(TimeDeathTest, MonthEndIncrementAborts) {
    const Calendar cal = Calendar::Gregorian;
    const Time jan31 = Date{2024, 1, 31}.to_time(cal);
    EXPECT_DEATH(jan31.add_months(cal, 1), "invalid date 2024-02-31");
    // Increments do not decompose: the same offset in one call is fine.
    EXPECT_EQ(jan31.add_months(cal, 2).to_date(cal), (Date{2024, 3, 31}));
    // add_years keeps month and day, so it aborts the same way off a leap day.
    const Time leap = Date{2024, 2, 29}.to_time(cal);
    EXPECT_DEATH(leap.add_years(cal, 1), "invalid date 2025-02-29");
}

// NO_CALENDAR is pure time spans: it has no dates, so both conversions abort
// rather than silently picking a calendar.
TEST(TimeDeathTest, NoCalendarHasNoDates) {
    EXPECT_DEATH(Time{}.to_date(Calendar::NoCalendar),
                 "TIM::Time::to_date: undefined for NO_CALENDAR");
    EXPECT_DEATH(Date{}.to_time(Calendar::NoCalendar),
                 "TIM::Date::to_time: undefined for NO_CALENDAR");
    EXPECT_DEATH(Time{}.add_months(Calendar::NoCalendar, 1),
                 "TIM::Time::add_months: undefined for NO_CALENDAR");
}

// A point before the calendar base is a legitimate span endpoint, but has no date.
TEST(TimeDeathTest, TimeBeforeBaseDateHasNoDate) {
    const Time before = Time{} - 1s;
    EXPECT_LT(before, Time{});
    EXPECT_EQ(Time{} - before, 1s);
    EXPECT_DEATH(before.to_date(Calendar::Gregorian),
                 "precedes the calendar base date \\(0001-01-01\\) by 1 s");
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
