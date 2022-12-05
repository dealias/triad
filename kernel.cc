/*                          T R I A D
An object-oriented C++ package for integrating initial value problems.
Copyright (C) 2000-2016 John C. Bowman (bowman@math.ualberta.ca)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <omp.h>

#include "options.h"
#include "kernel.h"
#include "fpu.h"

#include <sys/stat.h> // On sun machines this must come after xstream.h

using namespace Array;
using namespace std;

const char PROGRAM[]="TRIAD";
const char VERSION[]="1.46";
const int RestartVersion=4;
Real RestartProblemVersion=0;

// Global variables
size_t nout;
double t;
double last_dump=-1.0;
size_t iteration=0;
size_t invert_cnt=0;

VocabularyBase *Vocabulary;
ProblemBase *Problem;
IntegratorBase *Integrator;

// Kernel variables
static ::vector y;
static size_t ny;
static bool explicit_dt=0;
static bool testing=0;
static double cpu[ncputime],cpu0[ncputime],cpu_restart[ncputime];
static size_t final_iteration=0;
static size_t total_invert_cnt=0;

static const char *pname,*rname,*ptemp,*rtemp,*lname;
static ifstream fparam,fin;
static ofstream gparam,fdump,fstats,flock;

// Global vocabulary declarations and default values
size_t itmax=0;
size_t checkpoint=0;
double tmax=0.0;
double dt=0.0;
double tolmax=0.005;
double tolmin=0.003;
double stepfactor=2.0;
double stepnoninvert=stepfactor;
double dtmin=0.0;
double dtmax=0.0;
double tprecision=0.0;
double polltime=0.0;
size_t dynamic=1;
size_t digits=REAL_DIG;
size_t restart=0;
size_t output=0;
size_t hybrid=0;
size_t override=0;
size_t verbose=1;
size_t threads=0;

// Local vocabulary declarations and default values
size_t microsteps=1;
Real microfactor=1.0;
static double sample=0.0;
static size_t initialize=0;
static size_t clobber=0;
static const char *tmpdir=tempdir();

VocabularyBase::VocabularyBase()
{
  Vocabulary=this;

  VOCAB(itmax,(size_t) 0,SIZE_MAX,"Maximum number of iterations");
  VOCAB(microsteps,(size_t) 1,SIZE_MAX,"Number of iterations per output step");
  VOCAB(microfactor,0.0,REAL_MAX,"Microstep growth factor per output step");
  VOCAB(tmax,0.0,DBL_STD_MAX,"");
  VOCAB(dt,0.0,DBL_STD_MAX,"");
  VOCAB(dynamic,0,1,"");
  VOCAB(tolmax,0.0,DBL_STD_MAX,"");
  VOCAB(tolmin,0.0,DBL_STD_MAX,"");
  VOCAB(stepfactor,1.0,DBL_STD_MAX,"");
  VOCAB(stepnoninvert,1.0,DBL_STD_MAX,"");
  VOCAB(dtmin,0.0,DBL_STD_MAX,"");
  VOCAB(dtmax,0.0,DBL_STD_MAX,"");
  VOCAB(tprecision,0.0,DBL_STD_MAX,"");
  VOCAB(sample,0.0,DBL_STD_MAX,"Number of time units per output step");
  VOCAB(polltime,0.0,DBL_STD_MAX,"");
  VOCAB(hybrid,0,1,"");
  VOCAB(digits,1,SIZE_MAX,"");
  VOCAB_NODUMP(restart,0,1,"");
  VOCAB_NODUMP(initialize,0,1,"");
  VOCAB_NODUMP(clobber,0,1,"");
  VOCAB_NODUMP(override,0,1,"");
  VOCAB_NODUMP(verbose,0,4,"");
  VOCAB_NODUMP(threads,0,SIZE_MAX,"");
  VOCAB_NODUMP(run,null,null,"");
  VOCAB(checkpoint,(size_t) 0,SIZE_MAX,"");
  VOCAB(output,0,1,"");
  VOCAB_NOLIMIT(method,"");
  VOCAB_NOLIMIT(integrator,"");

  ProblemTable=new Table<ProblemBase>("method");
  IntegratorTable=new Table<IntegratorBase>("integrator");

  INTEGRATOR(Exact);
  INTEGRATOR(Euler);
  INTEGRATOR(Midpoint);
  INTEGRATOR(AB2);
  INTEGRATOR(ABM3);
  INTEGRATOR(PC);
  INTEGRATOR(LeapFrog);
  INTEGRATOR(RK1);
  INTEGRATOR(RK2);
  INTEGRATOR(RKPC);
  INTEGRATOR(RK3);
  INTEGRATOR(RK3C);
  INTEGRATOR(RK4);
  INTEGRATOR(RK5);
  INTEGRATOR(RK5DP);
  INTEGRATOR(SYM1);
  INTEGRATOR(SYM2);
}

void adjust_parameters(double& dt, double& dtmax, double& tmax,
		       size_t& itmax)
{
  if(dt == 0.0) {
    if(tmax == 0.0) tmax=1.0;
    dt=fabs(tmax/(itmax ? itmax : 1));
  }

  if(tmax == 0.0) tmax=DBL_STD_MAX;

  if(dtmax == 0.0) dtmax=DBL_STD_MAX;

  if(dt < dtmin) dt=dtmin;
  if(dt > dtmax) dt=dtmax;
}

int main(int argc, char *argv[])
{
  inform=mailuser;
  cout.precision(REAL_DIG);

  cout << newl << PROGRAM << " version " << VERSION
       << " [(C) John C. Bowman 2000]" << newl;

  cout << newl << "MACHINE: " << machine() << " [" << date() << "]" << newl;

  cout << newl << "PROBLEM: " << Vocabulary->Name() << " (version "
       << ProblemVersion << ")" << newl;

  cout << newl << "COMMAND LINE: ";
  for(int i=1; i < argc; i++) cout << argv[i] << " ";
  cout << endl;

  fpu_trap();

  Vocabulary->Sort();

  // Do a preliminary parse of command line to obtain parameter file name.
  for(int i=1; i < argc; i++) Vocabulary->Assign(argv[i],0);
  msg_override=override;

  // Allow time step to be overridden from command line (even on restarts).
  if(dt) explicit_dt=1;

  if(*run == 0) {run="test"; testing=1;}

  lname=Vocabulary->FileName(dirsep,"LOCK");
  rname=Vocabulary->FileName(dirsep,"restart");
  rtemp=Vocabulary->FileName(dirsep,"restart=");
  pname=Vocabulary->FileName(dirsep,"p");
  ptemp=Vocabulary->FileName(dirsep,"p=");

  fparam.open(pname);

  if(fparam) {
    string s;
    ostringstream buf;
    while(getline(fparam,s)) {
      buf << s << " ";
    }
    const string& S=buf.str();
    Vocabulary->Parse((char *) S.c_str());
    fparam.close();
  } else {
    if(!testing)
      msg(OVERRIDE_GLOBAL,"Parameter file %s could not be opened",
	  pname);
    errno=0;
    mkdir(Vocabulary->FileName(dirsep),0xFFFF);
  }

  for(int i=1; i < argc; i++) Vocabulary->Assign(argv[i]);
  adjust_parameters(dt,dtmax,tmax,itmax);
  cout << newl << "PARAMETERS:" << newl << newl;
  Vocabulary->List(cout);

  fpu_trap(!dynamic);
  cout.precision(digits);

  Problem=Vocabulary->NewProblem(method);
  Integrator=Vocabulary->NewIntegrator(integrator);

  if(!(restart || initialize)) {
    fin.open(rname);
    if(fin && !clobber && !testing)
      msg(OVERRIDE_GLOBAL,"Restart file %s already exists",rname);
    fin.close();
    errno=0;
  }

  if(restart) testlock();

  if(!restart && initialize) {
    const char *sname=Vocabulary->FileName(dirsep,"stat");
    fin.open(sname);
    if(fin && !clobber && !testing)
      msg(OVERRIDE_GLOBAL,"Statistics file %s already exists",sname);
    fin.close();
    errno=0;
  }

  if(checkpoint && tmpdir) {
    ostringstream buf;
    buf << tmpdir << dirsep << Vocabulary->FileName("","");
    const string& s=buf.str();
    mkdir(s.c_str(),0xFFFF);
  }

  SaveParameters();
  if(threads == 0) threads=omp_get_max_threads();
  cout << newl << "THREADS: " << threads << endl;
  t=0.0; nout=0;
  Problem->InitialConditions();

  Dimension(y,Problem->yVector());
  ny=Problem->Size();

  if(restart || initialize) read_init();
  if(!restart) Problem->Initialize();

  for(size_t i=0; i < 2; ++i) {
    open_output(gparam,(i == 0 ? "" : dirsep),"param.asy",0);
    Vocabulary->GraphicsDump(gparam,i);
    gparam.close();
  }

  Problem->Setup();
  Integrator->SetParam(tolmax,tolmin,stepfactor,stepnoninvert,dtmin,dtmax,
		       itmax,microsteps,microfactor,verbose,dynamic);

  Integrator->Allocator(*Problem,Problem->align);

  cout << newl << "INTEGRATING:" << endl;
  set_timer();

  Integrator->Integrate(t,tmax,dt,sample,iteration,nout);

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
  size_t i,ny0,nout0;
  double t0,dt0;
  const char *type=restart ? "restart" : "initialization";
  const char *ny_msg=
    "Current value of ny (%d) disagrees with value (%d) in file\n%s.";
  size_t init_version=0;

  ixstream finit(rname);
  if(finit) {
    cout << newl << "READING " << upcase(type) << " DATA FROM FILE " << rname;
    finit >> init_version >> nout0;
    if(init_version > 1) {
      finit >> RestartProblemVersion;
      cout << " (version " << RestartProblemVersion << ")";
    }
    cout << "." << endl;
    finit >> t0 >> dt0;

    if(init_version < 3) {
      size_t final_iteration0;
      finit >> final_iteration0;
      final_iteration=final_iteration0;
    } else finit >> final_iteration;

    for(i=0; i < ncputime; i++) finit >> cpu_restart[i];
    finit >> ny0;
    if(ny0 != ny) msg(OVERRIDE_GLOBAL,ny_msg,ny,ny0,rname);
    for(i=0; i < min(ny,ny0); i++) finit >> y[i];
    if(!finit)	msg(ERROR,"Cannot read from %s file %s",type,rname);
  } else msg(ERROR,"Could not open %s file %s",type,rname);

  if(!explicit_dt) dt=dt0;
  if(restart) {t=t0; nout=nout0; last_dump=t;}
}

void lock()
{
  if(testing) return;
  flock.open(lname);
  if(!flock) msg(WARNING,"Could not create lock file %s",lname);
}

void unlock()
{
  if(testing) return;
  if(flock) {
    flock.close();
    if(remove(lname) == -1)
      msg(WARNING,"Could not remove lock file %s",lname);
  }
}

void testlock()
{
  ifstream ftest;
  ftest.open(lname);
  if(ftest) {
    msg(OVERRIDE_GLOBAL,
	"Lock file %s exists.\nCheck output files",lname);
    if(remove(lname) == -1)
      msg(WARNING,"Could not remove lock file %s",lname);
  } else errno=0;
}

void dump(double t, size_t it, size_t final, double tmax)
{
  if((!restart || it > 0) && (tmax-t >= tprecision*tmax || final) &&
     t > last_dump) {
    lock();
    Problem->Output(it); last_dump=t;
    Problem->PrepareOutput(false);
    unlock();
  }
}

void SaveParameters() {
  fdump.open(ptemp);
  fdump.precision(REAL_DIG);
  Vocabulary->Dump(fdump);
  fdump.close();
  if(fdump) {
    if(rename(ptemp,pname))
      msg(WARNING,"Cannot rename parameter file %s",ptemp);
  }
  else msg(WARNING,"Cannot write to parameter file %s",ptemp);
}

static const size_t w=10;
static size_t e;

void set_timer()
{
  e=digits+5;
  open_output(fstats,dirsep,"stat");
  cputime(cpu0);
  if(!restart) {
    for(size_t i=0; i < ncputime; i++) cpu_restart[i]=0.0;
    fstats << "#" << setw(w) << "it"
	   << " " << setw(w) << "iteration"
	   << " " << setw(e) << "t"
	   << " " << setw(e) << "dt"
	   << " " << setw(w) << "invert_cnt"
	   << " " << setw(e) << "CPU"
	   << " " << setw(e) << "CHILD"
	   << " " << setw(e) << "SYS"
	   << " " << setw(w) << "memory"
	   << " " << setw(e) << "cpu" << endl;
  }
}

void statistics(double t, double dt, size_t it)
{
  if(restart && it == 0) return;
  static double lastcputime=0.0;

  static size_t last_iter=0;
  size_t iter=final_iteration+iteration;

  cputime(cpu);
  for(size_t i=0; i < ncputime; i++) cpu[i] -= cpu0[i];

  oxstream frestart(rtemp);
  if(frestart) {
    frestart << RestartVersion << newl << nout << newl
	     << ProblemVersion << newl
	     << t << newl << dt << newl << iter << newl;
    for(size_t i=0; i < ncputime; i++) frestart << cpu_restart[i]+cpu[i] << newl;
    frestart << ny << newl;
    for(size_t i=0; i < ny; i++) frestart << y[i] << newl;
    frestart.close();
    if(frestart) {
      if(checkpoint && it > 0 && (it-1) % checkpoint == 0) {
	ostringstream rcheck;
	if(tmpdir) rcheck << tmpdir << dirsep;
	rcheck << rname << "." << last_iter;
        const string& s=rcheck.str();
	if(tmpdir ? copy(rname,s.c_str()) : rename(rname,s.c_str()))
	  msg(WARNING,"Cannot copy %s to checkpoint file %s",
	      rname,s.c_str());
      }
      last_iter=iter;
      if(rename(rtemp,rname))
	msg(ERROR,"Cannot rename restart file %s",rtemp);
    } else {
      errno=0;
      frestart.open(rtemp);
      frestart.close();
      msg(SLEEP,"Cannot write to restart file %s",rtemp);
      statistics(t,dt,it); // Try again;
      return;
    }
  } else if(it == 0) msg(ERROR,"Could not open restart file %s",rtemp);

  lock();
  fstats << setw(w) << it << " " << iter << " " << setw(e) << t << " "
	 << setw(e) << dt << " " << setw(w) << invert_cnt << " ";

  total_invert_cnt += invert_cnt;
  invert_cnt=0;

  for(size_t i=0; i < ncputime; i++) {
    fstats << setw(e) << cpu_restart[i]+cpu[i] << " ";
  }
  fstats << setw(w) << memory() << " "
	 << setw(e) << cpu[0]-lastcputime << endl;
  lastcputime=cpu[0];

  if(!fstats) msg(WARNING,"Cannot write to statistics file");
  unlock();
}

const char *VocabularyBase::FileName(const char* delimiter, const char *suffix)
{
  ostringstream buf;
  buf << Directory();
  if(*delimiter) buf << run << delimiter;
  buf << suffix;
  const string& s=buf.str();
  return strdup(s.c_str());
}

void check_compatibility(const bool debug)
{
  const char *incompatible="Compiled with incompatible debugging modes";
  if(debug != DEBUG) msg(ERROR,incompatible);
}
