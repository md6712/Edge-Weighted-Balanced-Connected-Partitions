#ifndef Timer_H
#define Timer_H

#include <time.h>

struct Timer {

	
	long starttime;
	long endtime;
	long currenttime;

public:
	void setStartTime();
	void setEndTime();
	int calcElaspedTime();
	double calcElaspedTime_sec();
	int calcElaspedTime_curr();
	double calcElaspedTime_sec_curr();
};



#endif