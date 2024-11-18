#include "timer.h"

void Timer::setStartTime() {
	currenttime = starttime = clock();
}
void Timer::setEndTime() {
	currenttime = endtime = clock();
}
int Timer::calcElaspedTime() {
	return endtime - starttime;
}
double Timer::calcElaspedTime_sec() {
	return (double)(endtime - starttime) / CLK_TCK;
}

int Timer::calcElaspedTime_curr() {
	return currenttime - starttime;
}
double Timer::calcElaspedTime_sec_curr() {
	return (double)(currenttime - starttime) / CLK_TCK;
}