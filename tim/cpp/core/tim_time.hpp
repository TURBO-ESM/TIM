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
#include <cstdint>
#include <string>

namespace TIM {

/// @brief An exact, calendar-free span: the difference between two Time
/// points, and the increment for every unit of fixed length (seconds to days).
using Duration = std::chrono::duration<std::int64_t>;

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
/// @note Aborts if @p fms_type is outside 0..4.
Calendar calendar_from_fms(int fms_type);

/// @brief The length of a month in days.
/// @param cal The calendar to measure in.
/// @param year Calendar year.
/// @param month Month of year.
/// @return The number of days in that month.
/// @note Aborts if @p month is outside 1..12.
int days_in_month(Calendar cal, int year, int month);

class Time;

/// @brief A calendar date and time-of-day aggregate.
struct Date {
    int year = 1;    ///< Calendar year (>= 1; there is no year 0).
    int month = 1;   ///< Month of year, 1..12.
    int day = 1;     ///< Day of month, 1-based.
    int hour = 0;    ///< Hour of day, 0..23.
    int minute = 0;  ///< Minute of hour, 0..59.
    int second = 0;  ///< Second of minute, 0..59.

    /// Chronological ordering
    friend constexpr auto operator<=>(const Date&, const Date&) = default;

    /// @brief Render this date for diagnostics and error messages.
    /// @return The date as "YYYY-MM-DD hh:mm:ss".
    std::string to_string() const;

    /// @brief This date as a point on the time axis (FMS set_date).
    /// @param cal The calendar to read this date in.
    /// @return The resulting time.
    /// @note Aborts if this date does not exist in @p cal.
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

    /// @brief The exact span from the calendar base date to this point.
    /// @return That span, in seconds; negative if this point precedes the base date.
    constexpr Duration since_base() const noexcept { return since_base_; }

    /// @brief Chronological ordering.
    friend constexpr auto operator<=>(const Time&, const Time&) = default;

    // ---- calendar-free arithmetic ----
    //
    // A Duration is an exact span, so shifting a Time by one needs no calendar:
    // the result is the same point in every calendar. Only the conversions
    // below (to_date, add_months, add_years) have to know what a month is.
    //
    // Time deliberately has no operator* or operator/, which FMS provides on
    // time_type: scaling and dividing are span operations, and std::chrono
    // already supplies them on Duration.

    /// @brief Advance a point by a span.
    /// @param t The starting point.
    /// @param d The span to advance by; negative moves backwards.
    /// @return The point @p d after @p t.
    friend constexpr Time operator+(Time t, Duration d) noexcept { return Time{t.since_base_ + d}; }

    /// @brief Advance a point by a span (commuted form).
    /// @param d The span to advance by; negative moves backwards.
    /// @param t The starting point.
    /// @return The point @p d after @p t.
    friend constexpr Time operator+(Duration d, Time t) noexcept { return t + d; }

    /// @brief Move a point back by a span.
    /// @param t The starting point.
    /// @param d The span to retreat by; negative moves forwards.
    /// @return The point @p d before @p t.
    friend constexpr Time operator-(Time t, Duration d) noexcept { return Time{t.since_base_ - d}; }

    /// @brief The span between two points. Unlike FMS, which clamps the
    /// difference at zero, this is signed and needs no caller-side workaround.
    /// @param to The later point (by convention).
    /// @param from The earlier point (by convention).
    /// @return The exact span from @p from to @p to; negative if @p to precedes it.
    friend constexpr Duration operator-(Time to, Time from) noexcept {
        return to.since_base_ - from.since_base_;
    }

    /// @brief Advance this point in place by a span.
    /// @param d The span to advance by; negative moves backwards.
    /// @return This point, advanced.
    constexpr Time& operator+=(Duration d) noexcept { since_base_ += d; return *this; }

    /// @brief Move this point back in place by a span.
    /// @param d The span to retreat by; negative moves forwards.
    /// @return This point, moved back.
    constexpr Time& operator-=(Duration d) noexcept { since_base_ -= d; return *this; }

    // ---- calendar-dependent (FMS get_date / increment_date) ----

    /// @brief Convert this time to a date.
    /// @param cal The calendar to express the date in.
    /// @return The resulting date.
    /// @note Aborts if this point has no date in @p cal.
    Date to_date(Calendar cal) const;

    /// @brief This time advanced by @p n months, keeping the day of month.
    /// @param cal The calendar deciding what a month is.
    /// @param n Month count; negative allowed.
    /// @return The advanced point.
    /// @note The day of month is kept as-is, so this aborts whenever that day
    /// does not exist in the target month: Jan 31 + 1 month is 2024-02-31,
    /// which is fatal.
    /// @note Increments do not decompose: add_months(cal, 2) from Jan 31
    /// succeeds where add_months(cal, 1) twice aborts on the intermediate.
    /// Pass the full offset in one call rather than stepping.
    Time add_months(Calendar cal, int n) const;

    /// @brief This time advanced by @p n years, keeping the month and day.
    /// @param cal The calendar deciding what a year is.
    /// @param n Year count; negative allowed.
    /// @return The advanced point.
    /// @note As with add_months, the day of month is kept as-is: Feb 29 + 1
    /// year is 2025-02-29, which aborts.
    Time add_years(Calendar cal, int n) const;

private:
    Duration since_base_{};  ///< Signed exact span from the calendar base date (seconds).
};

}  // namespace TIM
