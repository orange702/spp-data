#include "gnsstime.h"
#include <cmath>

// ������ʱ�䣨�ꡢ�¡��ա�ʱ���֡��룩ת��ΪGPS�ܺ�������
// ����һ����³����ʵ�֣�ʹ���������������⸡�㾫������
void ymdhmsToGPS(int year, int month, int day, int hour, int min, double sec, int& week, double& sow) {
    // ���㷨���������գ�Julian Day���ļ���
    // �ο�: "Explanatory Supplement to the Astronomical Almanac"

    // ��������ת��Ϊ�������� (JD)
    // 1. ����·���1�»�2�£���Ϊ��һ���13�»�14��
    if (month <= 2) {
        month += 12;
        year--;
    }

    // 2. ����������
    // �����ʽ�Ǳ�׼�ģ����Ծ�ȷ����ӹ�Ԫǰ4713�꿪ʼ������
    double jd = floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day - 1524.5;
    jd += (hour / 24.0) + (min / 1440.0) + (sec / 86400.0);

    // GPSʱ�������� 1980��1��6��00:00:00 UTC����Ӧ���������� 2444244.5
    if (jd < 2444244.5) {
        // ���ʱ������GPS��ʼʱ�䣬�򷵻���Чֵ
        week = 0;
        sow = 0.0;
        return;
    }

    // 3. ����GPS��
    // (JD - 2444244.5) �Ǵ�GPS��ʼʱ�俪ʼ��������
    // ����7�õ�������
    week = floor((jd - 2444244.5) / 7.0);

    // 4. ����GPS������
    // (jd - 2444244.5) % 7.0 �ǵ�ǰ�ܵĵڼ��� (0=����, ..., 6=����)
    // ����һ������� (86400) �õ�������
    sow = fmod(jd - 2444244.5, 7.0) * 86400.0;
}


// ��GPS�ܺ�������ת��Ϊ����ʱ�䣨�˺����ڵ�ǰ��Ŀ�зǱ��裬�����������ԣ�
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