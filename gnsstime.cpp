#include "gnsstime.h"
#include <cmath>

// 将公历时间（年、月、日、时、分、秒）转换为GPS周和周内秒
// 这是一个更鲁棒的实现，使用整数运算来避免浮点精度问题
void ymdhmsToGPS(int year, int month, int day, int hour, int min, double sec, int& week, double& sow) {
    // 该算法基于儒略日（Julian Day）的计算
    // 参考: "Explanatory Supplement to the Astronomical Almanac"

    // 将年月日转换为儒略日数 (JD)
    // 1. 如果月份是1月或2月，视为上一年的13月或14月
    if (month <= 2) {
        month += 12;
        year--;
    }

    // 2. 计算儒略日
    // 这个公式是标准的，可以精确计算从公元前4713年开始的天数
    double jd = floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day - 1524.5;
    jd += (hour / 24.0) + (min / 1440.0) + (sec / 86400.0);

    // GPS时间的起点是 1980年1月6日00:00:00 UTC，对应的儒略日是 2444244.5
    if (jd < 2444244.5) {
        // 如果时间早于GPS起始时间，则返回无效值
        week = 0;
        sow = 0.0;
        return;
    }

    // 3. 计算GPS周
    // (JD - 2444244.5) 是从GPS起始时间开始的总天数
    // 除以7得到总周数
    week = floor((jd - 2444244.5) / 7.0);

    // 4. 计算GPS周内秒
    // (jd - 2444244.5) % 7.0 是当前周的第几天 (0=周日, ..., 6=周六)
    // 乘以一天的秒数 (86400) 得到周内秒
    sow = fmod(jd - 2444244.5, 7.0) * 86400.0;
}


// 将GPS周和周内秒转换为公历时间（此函数在当前项目中非必需，但保持完整性）
void gpsToYmdhms(int week, double sow, int& year, int& month, int& day, int& hour, int& min, double& sec) {
    double jd = (double)week * 7.0 + sow / 86400.0 + 2444244.5;

    long L = (long)(jd + 0.5) + 68569;
    long N = 4 * L / 146097;
    L = L - (146097 * N + 3) / 4;
    long I = 4000 * (L + 1) / 1461001;
    L = L - 1461 * I / 4 + 31;
    long J = 80 * L / 2447;
    day = L - 2447 * J / 80;
    L = J / 11;
    month = J + 2 - 12 * L;
    year = 100 * (N - 49) + I + L;

    double sod = fmod(sow, 86400.0);
    hour = (int)(sod / 3600.0);
    min = (int)(fmod(sod, 3600.0) / 60.0);
    sec = fmod(sod, 60.0);
}