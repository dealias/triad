/*                          T R I A D    
An object-oriented C++ package for integrating initial value problems.
Copyright (C) 1997 John C. Bowman (bowman@ipp-garching.mpg.de)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include "options.h"
#include "kernel.h"

#include <iomanip.h>

const char PROGRAM[]="TRIAD";
const char VERSION[]="1.1";

// Global variables
double t;
double last_dump=-1.0;
int iteration=0;
int invert_cnt=0;

VocabularyBase *Vocabulary;
ProblemBase *Problem;
IntegratorBase *Integrator;

// Kernel variables
static Var *y;
static int ny;
static int explicit_dt=0;
static int testing=0;
static double cpu[ncputime],cpu0[ncputime];
static final_iteration=0;
static int total_invert_cnt=0;

static char *pname,*rname,*ptemp,*rtemp,*lname;
static ifstream fparam,fin;
static ofstream fdump,fstat,flock;

// Global vocabulary declarations and default values
int itmax=-1;
double tmax=0.0;
double dt=0.0;
int Nmoment=0;
int dynamic=1;
double tolmax=0.005;
double tolmin=0.0;
double stepfactor=sqrt(2.0);
double stepnoninvert=stepfactor;
double dtmax=0.0;
int digits=REAL_DIG;
int restart=0;
double polltime=0.0;
int output=0;
int hybrid=0;
int override=0;
int verbose=1;
int checkpoint=0;

// Local vocabulary declarations and default values
static int microsteps=1;
static double sample=0.0;
static int initialize=0;
static int clobber=0;
static char *tmpdir=tempdir();

VocabularyBase::VocabularyBase()
{
	Vocabulary=this;
	
	VOCAB(itmax,0,INT_MAX);
	VOCAB(microsteps,1,INT_MAX);
	VOCAB(tmax,0.0,DBL_STD_MAX);
	VOCAB(dt,0.0,DBL_STD_MAX);
	VOCAB(dynamic,-INT_MAX,1);
	VOCAB(tolmax,0.0,DBL_STD_MAX);
	VOCAB(tolmin,0.0,DBL_STD_MAX);
	VOCAB(stepfactor,1.0,DBL_STD_MAX);
	VOCAB(stepnoninvert,1.0,DBL_STD_MAX);
	VOCAB(dtmax,0.0,DBL_STD_MAX);
	VOCAB(sample,0.0,DBL_STD_MAX);
	VOCAB(polltime,0.0,DBL_STD_MAX);
	VOCAB(hybrid,0,1);
	VOCAB(digits,1,INT_MAX);
	VOCAB_NODUMP(restart,0,1);
	VOCAB_NODUMP(initialize,0,1);
	VOCAB_NODUMP(clobber,0,1);
	VOCAB_NODUMP(override,0,1);
	VOCAB_NODUMP(verbose,0,4);
	VOCAB_NODUMP(run,"","");
	VOCAB(checkpoint,0,INT_MAX);
	VOCAB(output,0,1);
	VOCAB(method,"","");
	VOCAB(integrator,"","");
	
	ProblemTable=new Table<ProblemBase>("method");
	IntegratorTable=new Table<IntegratorBase>("integrator");
	
	INTEGRATOR(Exact);
	INTEGRATOR(Euler);
	INTEGRATOR(LeapFrog);
	INTEGRATOR(PC);
	INTEGRATOR(RK2);
	INTEGRATOR(RK4);
	INTEGRATOR(RK5);
}

void adjust_parameters(double& dt, double& dtmax, double& tmax, int& itmax)
{
	if(dt == 0.0) {
		if(tmax == 0.0) tmax=1.0;
		if(itmax == -1) itmax=100;
		dt=abs(tmax/(itmax ? itmax : 1));
	}
	
	if(tmax == 0.0) tmax=DBL_STD_MAX; 
	else if(itmax == -1) itmax=INT_MAX;

	if(dtmax == 0.0) dtmax=DBL_STD_MAX;
}	

int main(int argc, char *argv[])
{
	int i;
	cout.precision(REAL_DIG);
	
	cout << newl << PROGRAM << " version " << VERSION << 
		" [(C) John C. Bowman and B. A. Shadwick 1997]" << newl;
	
	cout << newl << "MACHINE: " << machine() << " [" << date() << "]" << newl;
	
	cout << newl << "PROBLEM: " << Vocabulary->Name() << newl;
	
	cout << newl << "COMMAND LINE: ";
	for(i=1; i < argc; i++) cout << argv[i] << " ";
	cout << endl;
	
	for(i=1; i < argc; i++) Vocabulary->Assign(argv[i]);
	
	// Allow time step to be overridden from command line (even on restarts).
	if(dt) explicit_dt=1; 
	
	if(*run == 0) {run="test"; testing=1; clobber=1;}
	
	lname=Vocabulary->FileName(dirsep,"LOCK");
	rname=Vocabulary->FileName(dirsep,"restart");
	rtemp=Vocabulary->FileName(dirsep,"restart=");
	pname=Vocabulary->FileName(dirsep,"p");
	ptemp=Vocabulary->FileName(dirsep,"p=");
	
	fparam.open(pname);
	
	if(fparam) {
		const int blocksize=80;
		char s[blocksize];
		while(!fparam.eof()) {
			strstream buf;
			while(1) {
				fparam.getline(s,blocksize);
				buf << s;
				if(fparam.eof() || !fparam.fail()) break;
				fparam.clear();
			}
			buf << ends;
			Vocabulary->Parse(buf.str());
		}
		fparam.close();
	} else {
		if(!testing) msg(ERROR,"Parameter file %s could not be opened",pname); 
		errno=0;
	}
	
	for(i=1; i < argc; i++) Vocabulary->Assign(argv[i]);
	adjust_parameters(dt,dtmax,tmax,itmax);
	cout << newl << "PARAMETERS:" << newl << newl;
	Vocabulary->List(cout);
	
	Problem=Vocabulary->NewProblem(method);
	Integrator=Vocabulary->NewIntegrator(integrator);
	
	if(!(restart || initialize)) {
		fin.open(rname);
		if(fin && !clobber) 
			msg(OVERRIDE,"Restart file %s already exists",rname);
		fin.close();	
		errno=0;
	}
	
	if(restart) testlock();

	if(!restart && initialize) {
		char *sname=Vocabulary->FileName(dirsep,"stat");
		fin.open(sname);
		if(fin && !clobber)
			msg(OVERRIDE,"Statistics file %s already exists",sname);
		fin.close();
		errno=0;
	}

	t=0.0;
	Problem->InitialConditions();
	y=Problem->Vector();
	ny=Problem->Size();
	if(restart || initialize) read_init();
	if(restart && dynamic < 0) dynamic=1;
	if(!restart) Problem->Initialize();
	
	if(!testing) {
		fdump.open(ptemp);
		Vocabulary->Dump(fdump);
		fdump.close();
		if(fdump) {
			if(rename(ptemp,pname)) 
				msg(WARNING,"Cannot rename parameter file %s",ptemp);
		}
		else msg(WARNING,"Cannot write to parameter file %s",ptemp);
	}
	
	Integrator->Allocate(ny);
	
	Integrator->SetParam(tolmax,tolmin,stepfactor,stepnoninvert,dtmax,itmax,
						 microsteps,verbose);
	
	cout << newl << "INTEGRATING:" << endl;
	set_timer();
	
	Integrator->Integrate(y,t,tmax,dt,sample);
	
	cputime(cpu);
	for(i=0; i < ncputime; i++) cpu[i] -= cpu0[i];

	Problem->FinalOutput();
	
	cout << newl;
	cout << "COMPUTATION TERMINATED after " << iteration <<	" iteration" <<
		PLURAL(iteration) << newl << "                       at t=" <<
		t << " with dt=" << dt << "." << newl; 
	if(total_invert_cnt)
		cout << "NONINVERTIBILE TRANSFORMATION encountered " <<
			total_invert_cnt << " time" << PLURAL(total_invert_cnt) <<
		"." << newl;
	
	cout << newl;
	cout << "INTEGRATION TIMING STATISTICS:" << newl <<	"            CPU = "
		 << cpu[0] << ", CHILD = " << cpu[1] <<  ", SYS = " << cpu[2] << newl;
	
	cout << newl << "DYNAMIC MEMORY ALLOCATED: " << memory() << " bytes ["
		 << date() << "]" << newl << endl; 
	
	if(!testing) mailuser("completed");
	return exit_signal;
}

void read_init()
{
	int i,ny0;
	double t0,dt0;
	char *type;
	if(restart) type="restart";
	else type="initialization";
	char *ny_msg=
		"Current value of ny (%d) disagrees with value (%d) in file\n%s";
	
	ixstream finit(rname);
	if(finit) {
		cout << newl << "READING " << upcase(type) << " DATA FROM FILE "
			 << rname << "." << endl; 
		finit >> t0 >> dt0 >> final_iteration;
		for(i=0; i < ncputime; i++) finit >> cpu[i];
		finit >> ny0;
		if(ny0 != ny) msg(OVERRIDE,ny_msg,ny,ny0,rname);
		for(i=0; i < min(ny,ny0); i++) finit >> y[i];
		if(!finit)	msg(ERROR,"Cannot read from %s file %s",type,rname);
	} else msg(ERROR,"Could not open %s file %s",type,rname);
	
	if(!explicit_dt) dt=dt0;
	if(restart) {t=t0; last_dump=t;}
}

void lock()
{
	flock.open(lname);
	if(!flock) {
		msg(WARNING,"Could not create lock file %s",lname);
		errno=0;
	} 
}

void unlock()
{
	if(flock) {
		flock.close();
		if(remove(lname) == -1) {
			msg(WARNING,"Could not remove lock file %s",lname);
			errno=0;
		}
	}
}

void testlock()
{
	ifstream ftest;
	ftest.open(lname);
	if(ftest) {
		msg(OVERRIDE,"Lock file %s exists.\nFiles may be corrupted",lname);
		flock.open(lname);
		unlock();
	} else errno=0;
}

void dump(int it, int final, double tmax) 
{
	if(!restart || it > 0) {
		if(!testing) lock();
		if((tmax-t >= 1.0E-6*tmax || final) && t > last_dump) {
			Problem->Output(it); last_dump=t;
		}
		statistics();
		if(!testing) unlock();
	}
	
	int iter=final_iteration+iteration;

	oxstream frestart(rtemp);
	if(frestart) {
		int i;
		frestart << t << newl << dt << newl << iter << newl;
		for(i=0; i < ncputime; i++) frestart << cpu[i] << newl;
		frestart << ny << newl;
		for(i=0; i < ny; i++) frestart << y[i] << newl;
		frestart.close();
		if(frestart) {
			if(checkpoint && it > 0 && (it-1) % checkpoint == 0) {
				strstream rcheck;
				if(tmpdir) rcheck << tmpdir << dirsep;
				rcheck << rname << "." << iter-microsteps << ends;
				if(rename(rname,rcheck.str()))
					msg(WARNING,"Cannot rename %s to checkpoint file %s",rname,
						rcheck.str());
			}
			if(rename(rtemp,rname))
				msg(WARNING,"Cannot rename restart file %s",rtemp);
		}
		else {
			errno=0;
			msg(WARNING,"Cannot write to restart file %s",rtemp);
		}		  
	} else if(it == 0) msg(ERROR,"Could not open restart file %s",rtemp);
}

static int w,e;

void set_timer()
{
	w=10;
	e=digits+5;
	open_output(fstat,dirsep,"stat");
	cputime(cpu0);
	if(restart) for(int i=0; i < ncputime; i++) cpu0[i] -= cpu[i];
	else fstat << setw(w) << "iteration" << " " << setw(e) << "t" << " " <<
		setw(e) << "dt" << " " << setw(w) << "invert_cnt" << " " <<
		setw(w) << "CPU" << " " << setw(w) << "CHILD" << " " <<
		setw(w) << "SYS" <<	endl;
}

void statistics()
{
	fstat << setw(w) << final_iteration+iteration << " " <<
		setw(e) << t << " " << setw(e) << dt << " " << setw(w) << 
		invert_cnt << " ";
	total_invert_cnt += invert_cnt;
	invert_cnt=0;
	cputime(cpu);
	for(int i=0; i < ncputime; i++) {
		cpu[i] -= cpu0[i];
		fstat << setw(w) << cpu[i] << " ";
	}
	fstat << endl;
}

char *VocabularyBase::FileName(const char* delimiter, const char *suffix)
{
	strstream buf;
	buf << Directory() << run << delimiter << suffix << ends;
	buf.rdbuf()->freeze();
	return buf.str();
}
