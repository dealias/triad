#include "kernel.h"

#include <iomanip.h>

const char PROGRAM[]="TRIAD";
const char VERSION[]="1.0";

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
static int override_dt=0;
static int testing=0;
static double cpu[ncputime],cpu0[ncputime];
static final_iteration=0;
static int total_invert_cnt=0;

static char *pname,*rname,*iname,*ptemp,*rtemp,*lname;
static ifstream fparam,fin;
static ofstream fdump,fstat,fout,flock;

// Global vocabulary declarations and default values
int itmax=0;
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
int discrete=0;

// Local vocabulary declarations and default values
static int microsteps=1;
static double sample=0.0;
static int initialize=0;
static int clobber=0;

static int average=0; // Obsolete
static char *approximation; // Obsolete

VocabularyBase::VocabularyBase()
{
	Vocabulary=this;
	
	VOCAB(itmax,0,INT_MAX);
	VOCAB(microsteps,1,INT_MAX);
	VOCAB(tmax,0.0,DBL_MAX);
	VOCAB(dt,0.0,DBL_MAX);
	VOCAB(dynamic,0,1);
	VOCAB(tolmax,0.0,DBL_MAX);
	VOCAB(tolmin,0.0,DBL_MAX);
	VOCAB(stepfactor,1.0,DBL_MAX);
	VOCAB(stepnoninvert,1.0,DBL_MAX);
	VOCAB(dtmax,0.0,DBL_MAX);
	VOCAB(sample,0.0,DBL_MAX);
	VOCAB(polltime,0.0,DBL_MAX);
	VOCAB(hybrid,0,1);
	VOCAB(digits,1,INT_MAX);
	VOCAB_NODUMP(restart,0,1);
	VOCAB_NODUMP(initialize,0,1);
	VOCAB_NODUMP(clobber,0,1);
	VOCAB_NODUMP(override,0,1);
	VOCAB_NODUMP(verbose,0,4);
	VOCAB_NODUMP(run,"","");
	VOCAB(output,0,1);
	VOCAB(discrete,0,1);
	VOCAB(problem,"","");
	VOCAB(integrator,"","");
	
	VOCAB_NODUMP(average,0,1); // Obsolete
	VOCAB_NODUMP(approximation,"",""); // Obsolete
		
	ProblemTable=new Table<ProblemBase>("problem",ProblemCompare,
										ProblemKeyCompare);
	IntegratorTable=new Table<IntegratorBase>("integrator",IntegratorCompare,
											  IntegratorKeyCompare);
	INTEGRATOR(Exact);
	INTEGRATOR(Euler);
	INTEGRATOR(PC);
	INTEGRATOR(RK2);
	INTEGRATOR(RK4);
	INTEGRATOR(RK5);
}

void adjust_parameters(double& dt, double& dtmax, double& tmax, int& itmax)
{
	if(dt == 0.0) {
		if(tmax == 0.0) tmax=1.0;
		if(itmax == 0) itmax=100;
		dt=abs(tmax/itmax);
	}
	
	if(tmax == 0.0) tmax=DBL_MAX; 
	else if(itmax == 0) itmax=INT_MAX;

	if(dtmax == 0.0) dtmax=DBL_MAX;
}	

int main(int argc, char *argv[])
{
	int i;

	cout.precision(REAL_DIG);
	
	cout << newl << PROGRAM << " version " << VERSION << 
		" [(C) John C. Bowman and B. A. Shadwick 1996]" << newl;
	
	cout << newl << "PROBLEM: " << Vocabulary->Name() << newl;
	
	cout << newl << "COMMAND LINE: ";
	for(i=1; i < argc; i++) cout << argv[i] << " ";
	cout << endl;
	
	for(i=1; i < argc; i++) Vocabulary->Assign(argv[i]);
	
	// Allow time step to be overridden from command line (even on restarts).
	if(dt) override_dt=1; 
	
	if(run == NULL) {run="test"; testing=1; clobber=1;}
	
	lname=Vocabulary->FileName(dirsep,"LOCK");
	rname=Vocabulary->FileName(dirsep,"restart");
	rtemp=Vocabulary->FileName(dirsep,"restart=");
	pname=Vocabulary->FileName(dirsep,"p");
	ptemp=Vocabulary->FileName(dirsep,"p=");
	
	fparam.open(pname);
	
	if(fparam) {
		const int size=256;
		char s[size],c;
		char *text="Parameter file %s contains a line longer than %d bytes";
		while(1) {
			fparam.get(s,size,'\n'); if(fparam.eof()) break;
			if(fparam.get(c) && c != '\n') msg(ERROR,text,pname,size);
			Vocabulary->Parse(s);
		}
		fparam.close();
	} else {
		if(!testing) msg(ERROR,"Parameter file %s could not be opened",pname); 
		errno=0;
	}
	
	for(i=1; i < argc; i++) Vocabulary->Assign(argv[i]);
	
	cout << newl << "PARAMETERS:" << newl << newl;
	Vocabulary->List(cout);
	
	
	if(!testing) {
		fdump.open(ptemp);
		Vocabulary->Dump(fdump);
		fdump.close();
		if(fdump.good()) rename(ptemp,pname);
		else msg(ERROR,"Cannot write to parameter file %s",ptemp);
	}
	
	Problem=Vocabulary->NewProblem(problem);
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
	}

	t=0.0;
	Problem->InitialConditions();
	y=Problem->Vector();
	ny=Problem->Size();
	if(restart || initialize) read_init();
	if(!restart) Problem->Initialize();
	
	Integrator->Allocate(ny);
	adjust_parameters(dt,dtmax,tmax,itmax);
	
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
	cout << endl;
	
	if(!testing) mailuser("completed");
	return exit_signal;
}

void read_init()
{
	ifstream finit;
	double t0,dt0;
	int ny0,formatted=0;
	char *ny_msg=
		"Current value of ny (%d) disagrees with value (%d) in file\n%s";
	
	iname=rname;
	finit.open(iname);
	errno=0;
	if(!finit) {
		formatted=1;
		iname=Vocabulary->FileName(dirsep,"restartf");
		finit.open(iname);
		if(!finit)
			msg(ERROR,"Initialization file %s could not be opened",iname);
	}
	
	cout << newl << "READING " << (restart ? "RESTART" : "INITIALIZATION") <<
		" DATA FROM FILE " << iname << "." << endl;

	if(formatted) {
		int i;
		finit >> t0 >> dt0;
		finit >> final_iteration;
		for(i=0; i < ncputime; i++) finit >> cpu[i];
		finit >> ny0;
		if(ny0 != ny) msg(OVERRIDE,ny_msg,ny,ny0,iname);
		for(i=0; i < min(ny,ny0); i++) finit >> y[i];
	} else {
		finit.read((char *) &t0,sizeof(double));
		finit.read((char *) &dt0,sizeof(double));
		finit.read((char *) &final_iteration,sizeof(int));
		finit.read((char *) cpu,sizeof(cpu));
		finit.read((char *) &ny0,sizeof(int));
		if(ny0 != ny) msg(OVERRIDE,ny_msg,ny,ny0,iname);
		finit.read((char *) y,min(ny,ny0)*sizeof(Var));
	}
	
	finit.close();
	if(!finit.good())
		msg(ERROR,"Cannot read from initialization file %s",iname);
	if(!override_dt) dt=dt0;
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
		if(remove(lname) == -1)
			msg(ERROR,"Could not remove lock file %s",lname);
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
	fdump.open(rtemp);
	if(fdump) {
		fdump.write((char *) &t,sizeof(double));
		fdump.write((char *) &dt,sizeof(double));
		fdump.write((char *) &iter,sizeof(int));
		fdump.write((char *) cpu,sizeof(cpu));
		fdump.write((char *) &ny,sizeof(int));
		fdump.write((char *) y,ny*sizeof(Var));
		fdump.close();
		if(fdump.good()) rename(rtemp,rname);
		else msg(WARNING,"Cannot write to restart file %s",rtemp);
	}
	else if(it == 0) msg(ERROR,"Dump file %s could not be opened",rtemp);

	if(output) {
		char *oname=Vocabulary->FileName(dirsep,"restartf");
		fout.open(oname);
		fout.precision(REAL_DIG);
		if(fout) {
			int i;
			fout << t << newl << dt << newl;
			fout << iter << newl;
			for(i=0; i < ncputime; i++) fout << cpu[i] << newl;
			fout << ny << newl;
			for(i=0; i < ny; i++) fout << y[i] << newl;
			fout.close();
		}
		if(!fout.good()) msg(WARNING,"Cannot write to output file %s",oname);
	}
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
	char *filename=new char[strlen(Abbrev())+strlen(run)+
	strlen(delimiter)+strlen(suffix)+2];
	sprintf(filename,"%s/%s%s%s",Abbrev(),run,delimiter,suffix);
	return filename;
}

void open_output(ofstream& fout, const char *delimiter, char *suffix,
				 int append)
{
	char *filename=Vocabulary->FileName(delimiter,suffix);
	if(append) fout.open(filename,ios::app); // Append to end of output file.
	else fout.open(filename);
	
	if(!fout) msg(ERROR,"Output file %s could not be opened",filename);
	fout.precision(digits);
}
