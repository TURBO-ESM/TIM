#pragma once
/**
 * @file tim_time.hpp
 * @brief TIM time concepts and calendar arithmetic.
 *
 * The three main concepts:
 *   - Duration: an exact, calendar-free span (a std::chrono duration).
 *   - Time:     a point on the time axis, measured from the calendar base.
 *   - Calendar: A tag used in converting Time to Date and vice versa.
 */

#include <chrono>
#include <compare>
#include <string>

namespace TIM {

/// @brief An exact, calendar-free span: the difference between two Time
/// points, and the increment for every unit of fixed length (seconds to days).
using Duration = std::chrono::seconds;

/// @brief The supported calendars. Enumerator values match the FMS
/// time_manager parameters, so a bind(C) bridge can pass FMS
/// get_calendar_type() straight to calendar_from_fms.
enum class Calendar : int {
    NoCalendar = 0,       ///< Pure time spans; no date arithmetic.
    ThirtyDayMonths = 1,  ///< Twelve 30-day months.
    Julian = 2,           ///< Every year divisible by 4 is leap.
    Gregorian = 3,        ///< Gregorian 4/100/400 leap rule.
    NoLeap = 4            ///< 365-day years, no leap days.
};

/// @brief Maps an FMS calendar-type integer to a Calendar.
/// @param fms_type FMS calendar parameter value (0..4).
/// @return The matching calendar.
Calendar calendar_from_fms(int fms_type);

/// @brief The length of a month in days.
/// @param cal The calendar to measure in.
/// @param year Calendar year.
/// @param month Month of year.
/// @return The number of days in that month.
int days_in_month(Calendar cal, int year, int month);

class Time;

/// @brief A calendar date and time-of-day agregate
struct Date {
    int year = 1;    ///< Calendar year (>= 1; there is no year 0).
    int month = 1;   ///< Month of year, 1..12.
    int day = 1;     ///< Day of month, 1-based.
    int hour = 0;    ///< Hour of day, 0..23.
    int minute = 0;  ///< Minute of hour, 0..59.
    int second = 0;  ///< Second of minute, 0..59.

    /// Chronological ordering
    friend constexpr auto operator<=>(const Date&, const Date&) = default;

    /// @brief "YYYY-MM-DD hh:mm:ss", for diagnostics and error messages.
    std::string to_string() const;

    /// @brief This date as a point on the time axis (FMS set_date).
    /// @param cal The calendar to read this date in.
    /// @return The resulting time.
    Time to_time(Calendar cal) const;
};

/// @brief A point on the time axis: the exact span from a calendar's base date
/// (the start of year 1) to that point. Time is calendar-free.
class Time {
public:
    /// @brief The calendar base date itself (a zero span).
    constexpr Time() = default;

    /// @brief The point @p since_base after the calendar base date.
    /// @param since_base The exact span from the calendar base date (in seconds)
    constexpr explicit Time(Duration since_base) noexcept : since_base_(since_base) {}

    /// @brief The exact span from the calendar base date to this point (in seconds).
    constexpr Duration since_base() const noexcept { return since_base_; }

    /// @brief Chronological ordering.
    friend constexpr auto operator<=>(const Time&, const Time&) = default;

    // ---- calendar-free arithmetic ----

    friend constexpr Time operator+(Time t, Duration d) noexcept { return Time{t.since_base_ + d}; }
    friend constexpr Time operator+(Duration d, Time t) noexcept { return t + d; }
    friend constexpr Time operator-(Time t, Duration d) noexcept { return Time{t.since_base_ - d}; }
    friend constexpr Duration operator-(Time to, Time from) noexcept {
        return to.since_base_ - from.since_base_;
    }
    constexpr Time& operator+=(Duration d) noexcept { since_base_ += d; return *this; }
    constexpr Time& operator-=(Duration d) noexcept { since_base_ -= d; return *this; }

    // ---- calendar-dependent (FMS get_date / increment_date) ----

    /// @brief Convert this time to a date.
    /// @param cal The calendar to express the date in.
    /// @return The resulting date.
    Date to_date(Calendar cal) const;

    /// @brief This time advanced by @p n months, keeping the day of month.
    /// @param cal The calendar deciding what a month is.
    /// @param n Month count; negative allowed.
    /// @return The advanced point.
    Time add_months(Calendar cal, int n) const;

    /// @brief This time advanced by @p n years, keeping the month and day.
    /// @param cal The calendar deciding what a year is.
    /// @param n Year count; negative allowed.
    /// @return The advanced point.
    Time add_years(Calendar cal, int n) const;

private:
    Duration since_base_{};  ///< Signed exact span from the calendar base date (seconds).
};

}  // namespace TIM
