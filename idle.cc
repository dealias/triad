#include <unistd.h>
#include "kernel.h"

// return CPU idle time in seconds

int idletime()
{	
	int idle,hh,mm,ss;
	FILE *idlefile;

	system("idletime > .idletime");
	if(idlefile=fopen(".idletime","r"))	{
		if(fscanf(idlefile,"%d:%d:%d\n",&hh,&mm,&ss) == 3)
			idle=(hh*60+mm)*60+ss;
		else idle=0;
		fclose(idlefile);
	}
	return idle;
}

// Check whether load average is too high. If so, sleep for a while.

void poll()
{
	double avg=1.0;
	int newidle=0,idle=0;
	const double sleeptime=polltime;
	char s[20];
	FILE *loadavg;
	
	system("uptime > .loadavg");
	if(loadavg=fopen(".loadavg","r")) {
		while(fscanf(loadavg,"%s",s) == 1) if(strcmp(s,"average:")==0) break;
		fscanf(loadavg,"%s",s);
		avg=atof(s);
		fclose(loadavg);
	}

	if(avg > 1.5) while(!idle || (newidle-idle) < 0.6*sleeptime) {
		idle=idletime();
		sleep(sleeptime);
		newidle=idletime();
	}
}


