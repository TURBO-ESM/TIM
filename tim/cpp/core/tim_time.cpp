/**
 * @file tim_time.cpp
 * @brief TIM time concepts and calendar arithmetic implementation.
 * Analogue of FMS time_manager.F90.
 */

#include "tim_time.hpp"

#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>

#include "tim_abort.hpp"

namespace TIM {

namespace {

constexpr int month_lengths[13] = {0,  31, 28, 31, 30, 31, 30,
                                   31, 31, 30, 31, 30, 31};
constexpr std::int64_t max_year = std::numeric_limits<int>::max();

constexpr std::int64_t seconds_per_day = 86400;
constexpr std::int64_t julian_era_days = 4 * 365 + 1;                       // 1461
constexpr std::int64_t gregorian_century_days = 100 * 365 + 24;             // 36524
constexpr std::int64_t gregorian_era_days = 4 * gregorian_century_days + 1;  // 146097

// The calendar's name, for diagnostics.
const char* calendar_name(const Calendar cal) {
    switch (cal) {
        case Calendar::NoCalendar:      return "NO_CALENDAR";
        case Calendar::ThirtyDayMonths: return "THIRTY_DAY_MONTHS";
        case Calendar::Julian:          return "JULIAN";
        case Calendar::Gregorian:       return "GREGORIAN";
        case Calendar::NoLeap:          return "NOLEAP";
    }
    return "UNKNOWN";
}

// Fortran's floor(a/real(b)) and modulo(a,b): the FMS increment_date month
// arithmetic, time_manager.F90:2551-2552. C's / and % truncate toward zero, so
// a negative month offset would collapse toward year 0 instead of borrowing a
// year.

constexpr std::int64_t floor_div(const std::int64_t a, const std::int64_t b) {
    const std::int64_t q = a / b;
    const std::int64_t r = a % b;
    return (r != 0 && ((r < 0) != (b < 0))) ? q - 1 : q;
}

constexpr std::int64_t floor_mod(const std::int64_t a, const std::int64_t b) {
    return a - floor_div(a, b) * b;
}

// Calendar-independent field ranges only (e.g., day 1..31)
constexpr bool valid_ranges(const Date& d) {
    return d.second >= 0 && d.second <= 59 && d.minute >= 0 && d.minute <= 59 &&
           d.hour >= 0 && d.hour <= 23 && d.day >= 1 && d.day <= 31 &&
           d.month >= 1 && d.month <= 12 && d.year >= 1;
}

// ---- calendar shape (year >= 1 everywhere) ----

// Only the two Julian-family calendars have leap years at all.
bool is_leap_year(const Calendar cal, const int year) {
    switch (cal) {
        case Calendar::Gregorian:
            return year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
        case Calendar::Julian: return year % 4 == 0;
        case Calendar::NoLeap:
        case Calendar::ThirtyDayMonths:
        case Calendar::NoCalendar: break;
    }
    return false;
}

// Serial day (0 = 0001-01-01) of the first day of year y: 365 (or 360) days
// per elapsed year plus the leap days among years [1, y-1].
std::int64_t days_before_year(const Calendar cal, const std::int64_t y) {
    const std::int64_t n = y - 1;
    switch (cal) {
        case Calendar::ThirtyDayMonths: return 360 * n;
        case Calendar::Gregorian:       return 365 * n + n / 4 - n / 100 + n / 400;
        case Calendar::Julian:          return 365 * n + n / 4;
        default:                        return 365 * n;  // NoLeap
    }
}

// Days from Jan 1 of `year` to the first of `month`. days_in_month is the one
// rule for month lengths, so there is no second table to keep in step with it.
std::int64_t days_before_month(const Calendar cal, const int year, const int month) {
    std::int64_t days = 0;
    for (int m = 1; m < month; ++m) days += days_in_month(cal, year, m);
    return days;
}

struct YearDoy {
    std::int64_t year;  // calendar year (>= 1)
    std::int64_t doy;   // 1-based day-of-year, 1..366
};

// Serial day -> (year, day-of-year), by the closed-form era algebra: divide the
// day-of-era by a common year, then correct so the extra days that leap years
// contribute never push the year-of-era past the end of its era.
YearDoy split_year(const Calendar cal, const std::int64_t days) {
    switch (cal) {
        case Calendar::ThirtyDayMonths:
            return {days / 360 + 1, days % 360 + 1};
        case Calendar::Julian: {
            const std::int64_t era = days / julian_era_days;
            const std::int64_t doe = days % julian_era_days;
            const std::int64_t yoe = (doe - doe / (julian_era_days - 1)) / 365;
            return {era * 4 + yoe + 1, doe - 365 * yoe + 1};
        }
        case Calendar::Gregorian: {
            const std::int64_t era = days / gregorian_era_days;
            const std::int64_t doe = days % gregorian_era_days;
            const std::int64_t yoe = (doe - doe / (4 * 365)
                                          + doe / gregorian_century_days
                                          - doe / (gregorian_era_days - 1)) / 365;
            return {era * 400 + yoe + 1, doe - (365 * yoe + yoe / 4 - yoe / 100) + 1};
        }
        default:  // NoLeap
            return {days / 365 + 1, days % 365 + 1};
    }
}


// Shift time t by whole calendar years and months, keeping day and clock time.
Time add_calendar_interval(const Time t, const Calendar cal, const std::int64_t years,
                           const std::int64_t months, const char* const who) {
    if (cal == Calendar::NoCalendar)
        TIM::abort(std::string("TIM::Time::") + who + ": undefined for NO_CALENDAR.");
    Date d = t.to_date(cal);
    const std::int64_t m = d.month + months;
    const std::int64_t y = d.year + floor_div(m - 1, 12) + years;
    if (y < 1 || y > max_year)
        TIM::abort(std::string("TIM::Time::") + who + ": the result falls outside the "
                  "representable years (1.." + std::to_string(max_year) + ").");
    d.year = static_cast<int>(y);
    d.month = static_cast<int>(floor_mod(m - 1, 12)) + 1;
    return d.to_time(cal);
}

}  // namespace

std::string Date::to_string() const {
    // todo: switch to std::format once it is available across all compilers
    char buf[40];
    std::snprintf(buf, sizeof buf, "%04d-%02d-%02d %02d:%02d:%02d", year, month, day,
                  hour, minute, second);
    return buf;
}

Calendar calendar_from_fms(const int fms_type) {
    if (fms_type < 0 || fms_type > 4)
        TIM::abort("TIM::calendar_from_fms: invalid FMS calendar type " +
                  std::to_string(fms_type) + "; expected 0..4.");
    return static_cast<Calendar>(fms_type);
}

int days_in_month(const Calendar cal, const int year, const int month) {
    if (month < 1 || month > 12)
        TIM::abort("TIM::days_in_month: month out of range: " +
                  std::to_string(month) + ".");
    if (cal == Calendar::ThirtyDayMonths) return 30;
    if (month == 2 && is_leap_year(cal, year)) return 29;
    return month_lengths[month];
}

Time Date::to_time(const Calendar cal) const {
    if (cal == Calendar::NoCalendar)
        TIM::abort("TIM::Date::to_time: undefined for NO_CALENDAR; a time span "
                  "under NO_CALENDAR has no date.");
    if (!valid_ranges(*this) || day > days_in_month(cal, year, month))
        TIM::abort("TIM::Date::to_time: invalid date " + to_string() +
                  " in calendar " + calendar_name(cal) + ".");
    const std::int64_t days = days_before_year(cal, year) +
                              days_before_month(cal, year, month) + day - 1;
    const std::int64_t seconds = second + 60 * (minute + 60 * hour);
    return Time{Duration{days * seconds_per_day + seconds}};
}

Date Time::to_date(const Calendar cal) const {
    if (cal == Calendar::NoCalendar)
        TIM::abort("TIM::Time::to_date: undefined for NO_CALENDAR; a time span "
                  "under NO_CALENDAR has no date.");
    const std::int64_t total = since_base_.count();
    if (total < 0)
        TIM::abort("TIM::Time::to_date: the time precedes the calendar base date "
                  "(0001-01-01) by " + std::to_string(-total) + " s; there is no year 0.");
    const std::int64_t days = total / seconds_per_day;
    const int seconds = static_cast<int>(total % seconds_per_day);  // 0..86399

    const auto [year, doy] = split_year(cal, days);
    if (year > max_year)
        TIM::abort("TIM::Time::to_date: the time falls beyond the representable "
                  "years (1.." + std::to_string(max_year) + ").");

    Date out;
    out.year = static_cast<int>(year);
    out.day = static_cast<int>(doy);
    for (out.month = 1; out.month <= 12; ++out.month) {
        const int month_days = days_in_month(cal, out.year, out.month);
        if (out.day <= month_days) break;
        out.day -= month_days;
    }
    if (out.month > 12)
        TIM::abort("TIM::Time::to_date: internal error: day-of-year " +
                   std::to_string(doy) + " is past the end of year " +
                   std::to_string(out.year) + " in calendar " + calendar_name(cal) + ".");
    out.hour = seconds / 3600;
    out.minute = seconds % 3600 / 60;
    out.second = seconds % 60;
    return out;
}

Time Time::add_months(const Calendar cal, const int n) const {
    return add_calendar_interval(*this, cal, 0, n, "add_months");
}

Time Time::add_years(const Calendar cal, const int n) const {
    return add_calendar_interval(*this, cal, n, 0, "add_years");
}

}  // namespace TIM
