#pragma once

// ������ʱ��ת��ΪGPS�ܺ�������
void ymdhmsToGPS(int year, int month, int day, int hour, int min, double sec, int& week, double& sow);

// ��GPS�ܺ�������ת��Ϊ����ʱ��
void gpsToYmdhms(int week, double sow, int& year, int& month, int& day, int& hour, int& min, double& sec);