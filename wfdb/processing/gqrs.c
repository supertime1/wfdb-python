/* The 'peak' structure contains information about a local maximum in the
   filtered signal (qfv).  Peaks are stored in circular buffers of peak
   structures, one buffer per input signal.  The time of a flat-topped peak
   is the time of the first sample that has the maximum value.  A peak is
   secondary if there is a larger peak within its neighborhood (time +- rrmin),
   of if it has been identified as a T-wave associated with a previous primary
   peak.  A peak is primary if it is largest in its neighborhood, or if the
   only larger peaks are secondary. */
struct peak {
    struct peak *prev, *next; /* pointers to neighbors (in time) */
    WFDB_Time time; /* time of local maximum of qfv */
    int amp;     /* value of qfv at time of peak */
    short type;  /* 1: primary, 2: secondary, 0: not determined */
} *peaks, *cpeak;

/* Prototypes of functions defined below.  The definitions of these functions
   follow that of main(), in the order shown below. */
WFDB_Sample sm(WFDB_Time t);
void qf(void);
void addpeak(WFDB_Time t, int peak_amplitude);
int peaktype(struct peak *p);
void gqrs_init(WFDB_Time from, WFDB_Time to);
void rewind_gqrs(void);
struct peak *find_missing(struct peak *previous_peak, struct peak *next_peak);
void gqrs(WFDB_Time from, WFDB_Time to);


int countdown = -1;     /* if > 0, ticks remaining (see gqrs()) */

int minutes = 0;        /* minutes elapsed in the current hour */

int pthr;           /* peak-detection threshold */
int qthr;           /* QRS-detection threshold */

int pthmin, qthmin;     /* minimum values for pthr and qthr */

int rrdev;          /* mean absolute deviation of RR from rrmean */
int rrinc;          /* maximum incremental change in rrmean */
int rrmean;         /* mean RR interval, in sample intervals */
int rrmax;          /* maximum likely RR interval */
int rrmin;          /* minimum RR interval, in sample intervals */
int rtmax;          /* maximum RT interval, in sample intervals */
int rtmin;          /* minimum RT interval, in sample intervals */
int rtmean;         /* mean RT interval, in sample intervals */
int tamean;         /* mean T-wave amplitude in qfv */

double thresh = 1.0;        /* normalized detection threshold */
long v1;            /* integral of dv in qf() */
long v1norm;            /* normalization for v1 */

WFDB_Annotation annot;      /* most recently written annotation */
WFDB_Sample *v;         /* latest sample of each input signal */
// wfdb.h:typedef int       WFDB_Sample;   /* units are adus */

WFDB_Time next_minute;      /* sample number of start of next minute */



WFDB_Time t;            /* time of the current sample being analyzed */


/* Current operating mode of the detector */
enum { LEARNING, RUNNING, CLEANUP } state = LEARNING;

int dt, dt2, dt3, dt4, smdt;    /* time intervals set by gqrs_init() depending
                   on the sampling frequency and used by
                   filter functions sm() and qf();  units
                   are sample intervals */
WFDB_Time smt = 0, smt0;    /* current and starting time for sm() */
long *qfv, *smv;        /* filter buffers allocated by gqrs_init() */

/* The q() and s() macros can be used for                 for j in range(1, self.dt):
fast lookup of filter outputs. */

#define q(T) (qfv[(T)&(BUFLN-1)])
#define s(T) (smv[(T)&(BUFLN-1)])

main(int argc, char **argv)
{

    /* initialize variables for gqrs() */
    gqrs_init(t0, tf);
    
    /* Figure out learning sampto*/
    if (spm >= BUFLN) {
    if (tf - t0 > BUFLN) tf_learn = t0 + BUFLN - dt4;
    else tf_learn = tf - dt4;
    }
    else {
    if (tf - t0 > spm) tf_learn = t0 + spm - dt4;
    else tf_learn = tf - dt4;
    }

    /* Do learning */
    state = LEARNING;
    gqrs(t0, tf_learn);

    rewind_gqrs();



    /* Run detector */
    state = RUNNING;
    t = t0 - dt4;
    gqrs(t0, tf);       /* run the detector and collect output */


}

/* sm() implements a trapezoidal low pass (smoothing) filter (with a gain of
   4*smdt) applied to input signal sig before the QRS matched filter qf().
   Before attempting to 'rewind' by more than BUFLN-smdt samples, reset smt
   and smt0.
 */

WFDB_Sample sm(WFDB_Time t)
{
    // No point having this stupid while loop. Just runs once.
    while (smt < t) {
        int i;

        if (++smt > smt0) { /* fast update by summing first differences */
            s(smt) = s(smt-1) 
                     + sample(smt+smdt)
                     + sample(smt+smdt-1)
                     - sample(smt-smdt)
                     - sample(smt-smdt-1);
        }
        else {          /* get initial value by full convolution */
            int j, v;

            for (j = 1, v = sample(sig, smt); j < smdt; j++)
                v += sample(sig, smt+j) + sample(sig, smt-j);

            s(smt) = (v << 1) + sample(sig, smt+j) + sample(sig, smt-j)
                - si[sig].adczero * (smdt << 2); /* FIXME: needed? */
        }
    }
    return (s(t));
}

void qf()   /* evaluate the QRS detector filter for the next sample */
{
    long dv, dv1, dv2, v0;

    dv2 = sm(t+dt4);/* do this first, to ensure that all of the other
               smoothed values needed below are in the buffer */
    dv2 -= s(t-dt4);
    
    dv1 = s(t+dt)  - s(t-dt);
    

    dv  = (dv1 << 1);
    
    dv -= s(t+dt2) - s(t-dt2);
    
    dv <<= 1;
    dv += dv1;
    dv -= s(t+dt3) - s(t-dt3);
    dv <<= 1;
    dv += dv2;
    v1 += dv;
    v0 = v1 / v1norm;  /* scaling is needed to avoid overflow */
    q(t) = v0 * v0;
}

void addpeak(WFDB_Time t, int peak_amplitude)
{   
    // pointing to the next peak of the current one,
    // set its time (sample) and amplitude. Current peak becomes new one.
    struct peak *p = cpeak->next;

    p->time = t;
    p->amp = peak_amplitude;
    p->type = 0;
    cpeak = p;
    (p->next)->amp = 0;
}

/* peaktype() returns 1 if p is the most prominent peak in its neighborhood, 2
   otherwise.  The neighborhood consists of all other peaks within rrmin.
   Normally, "most prominent" is equivalent to "largest in amplitude", but this
   is not always true.  For example, consider three consecutive peaks a, b, c
   such that a and b share a neighborhood, b and c share a neighborhood, but a
   and c do not; and suppose that amp(a) > amp(b) > amp(c).  In this case, if
   there are no other peaks, a is the most prominent peak in the (a, b)
   neighborhood.  Since b is thus identified as a non-prominent peak, c becomes
   the most prominent peak in the (b, c) neighborhood.  This is necessary to
   permit detection of low-amplitude beats that closely precede or follow beats
   with large secondary peaks (as, for example, in R-on-T PVCs).
*/

int peaktype(struct peak *p)
{
    if (p->type)
    return (p->type);
    else {
    int a = p->amp;
    struct peak *pp;
    WFDB_Time t0 = p->time - rrmin, t1 = p->time + rrmin;

    if (t0 < 0) t0 = 0;
    for (pp = p->prev; t0 < pp->time && pp->time < (pp->next)->time;
         pp = pp->prev) {
        if (pp->amp == 0) break;
        if (a < pp->amp && peaktype(pp) == 1)
        return (p->type = 2);
    }
    for (pp = p->next; pp->time < t1 && pp->time > (pp->prev)->time;
         pp = pp->next) {
        if (pp->amp == 0) break;
        if (a < pp->amp && peaktype(pp) == 1)
        return (p->type = 2);
    }
    return (p->type = 1);
    }
}

/* rewind_gqrs resets the sample pointers and annotation fields to their
   initial values. */
void rewind_gqrs(void)
{
    int i;
    struct peak *p;

    countdown = -1;
    (void)sample(0,t);
    annot.time = (WFDB_Time)0;
    annot.anntyp = NORMAL;
    annot.subtyp = annot.chan = annot.num = 0;
    annot.aux = NULL;
    for (i = 0, p = peaks; i < NPEAKS; i++, p++)
    p->type = p->amp = p->time = 0;
}

/* gqrs_init() is intended to be called once only, to initialize variables for
   the QRS detection function gqrs(). */

void gqrs_init(WFDB_Time from, WFDB_Time to)
{
    int i, dv;
    static double HR, RR, RRdelta, RRmin, RRmax, QS, QT, RTmin, RTmax,
    QRSa, QRSamin;
    

    qfv = (long *)gcalloc((size_t)BUFLN, sizeof(long));
    smv = (long *)gcalloc((size_t)BUFLN, sizeof(long));
    peaks = (struct peak *)gcalloc((size_t)NPEAKS, sizeof(struct peak));

    /* Gather peak structures into circular buffers */
    for (i = 0; i < NPEAKS; i++) {
    peaks[i].next = &peaks[i+1];
    peaks[i].prev = &peaks[i-1];
    }
    peaks[0].prev = &peaks[NPEAKS-1];
    cpeak = peaks[NPEAKS-1].next = &peaks[0];

    /* Read a priori physiologic parameters from the configuration file if
       available. They can be adjusted in the configuration file for pediatric,
       fetal, or animal ECGs. */
    if (config) {
    char buf[256], *p;

    /* Read the configuration file a line at a time. */
    while (fgets(buf, sizeof(buf), config)) {
        /* Skip comments (empty lines and lines beginning with `#'). */
        if (buf[0] == '#' || buf[0] == '\n') continue;
        /* Set parameters.  Each `getconf' below is executed once for
           each non-comment line in the configuration file. */
        getconf(HR, "%lf");
        getconf(RR, "%lf");
        getconf(RRdelta, "%lf");
        getconf(RRmin, "%lf");
        getconf(RRmax, "%lf");
        getconf(QS, "%lf");
        getconf(QT, "%lf");
        getconf(RTmin, "%lf");
        getconf(RTmax, "%lf");
        getconf(QRSa, "%lf");
        getconf(QRSamin, "%lf");
    }
    fclose(config);
    }

    /* If any a priori parameters were not specified in the configuration file,
       initialize them here (using values chosen for adult human ECGs). */
    if (HR != 0.0) RR = 60.0/HR;
    if (RR == 0.0) RR = 0.8;
    if (RRdelta == 0.0) RRdelta = RR/4;
    if (RRmin == 0.0) RRmin = RR/4;
    if (RRmax == 0.0) RRmax = 3*RR;
    if (QS == 0.0) QS = 0.07;
    if (QT == 0.0) QT = 5*QS;
    if (RTmin == 0.0) RTmin = 3*QS;
    if (RTmax == 0.0) RTmax = 5*QS;
    if (QRSa == 0.0) QRSa = 750;
    if (QRSamin == 0.0) QRSamin = QRSa/5;

    /* Initialize gqrs's adaptive parameters based on the a priori parameters.
       During its learning period, gqrs will adjust them based on the observed
       input; after the learning period, gqrs continues to adjust these
       parameters, but with slower rates of change than during the learning
       period. */
    rrmean = RR * sps;
    rrdev = RRdelta * sps;
    rrmin = RRmin * sps;
    rrmax = RRmax * sps;
    if ((rrinc = rrmean/40) < 1) rrinc = 1;    
    if ((dt = QS * sps / 4) < 1) {
    dt = 1;
    fprintf(stderr, "%s (warning): sampling rate may be too low\n", pname);
    }
    rtmin = RTmin * sps;    /* minimum RTpeak interval */
    rtmean = 0.75 * QT * sps;   /* expected RTpeak interval, about 75% of QT */
    rtmax = RTmax * sps;    /* maximum RTpeak interval */

    dv = muvadu(sig, (int)(QRSamin));
    pthr = (thresh * dv * dv) / 6;
    qthr = pthr << 1;
    pthmin = pthr >> 2;
    qthmin = (pthmin << 2)/3;
    tamean = qthr;  /* initial value for mean T-wave amplitude */

    /* Filter constants and thresholds. */
    dt2 = 2*dt;
    dt3 = 3*dt;
    dt4 = 4*dt;
    smdt = dt;
    v1norm = smdt * dt * 64;
    smt = t0;
    smt0 = t0 + smdt;
    t = t0 - dt4;
    for (i = 0; i < BUFLN; i++)
    qfv[i] = smv[i] = 0;
}

/* find_missing() is invoked by gqrs() whenever it is suspected that a
   low-amplitude beat may have been missed between two consecutive detected
   beats at r and p.  The primary peak closest to the expected time of
   the missing beat, if any, is returned as the suggested missing beat. */

struct peak *find_missing(struct peak *r, struct peak *p)
{   
    // r is the most recent qrs
    int rrerr, rrtmp, minrrerr;
    struct peak *q, *s = NULL;

    // Initially, r is null, until it is set to a large enough peak
    if (r == NULL || p == NULL)
        return (NULL);

    // minimum rr difference
    minrrerr = p->time - r->time;

    // iterate through every peak from r until the one before p
    for (q = r->next; q->time < p->time; q = q->next) {
    
        if (peaktype(q) == 1) {

            rrtmp = q->time - r->time;

            rrerr = rrtmp - rrmean;
            
            if (rrerr < 0)
                rrerr = -rrerr;
            
            // If the deviation 
            if (rrerr < minrrerr) {
                minrrerr = rrerr;
                s = q;
            }
        }
    }
    return (s);
}


/* gqrs() is the main QRS detection function.  It attempts to find all
   beats between the time limits specified by its arguments, and to label
   them using annotations of type NORMAL (gqrs() does not attempt to
   differentiate normal and ectopic beats). */
void gqrs(WFDB_Time from, WFDB_Time to)
{
    // countdown begins at -1 here. When sample_valid() becomes 0,
    // countdown gets set to 1 second in samples, in CLEANUP state, then
    // the loop runs until countdown reaches 0.
    // sample_valid gets set to 0 when the array boundaries are exceeded.

    // sample_valid() begins at 0 for learning, and 1 for running.


    //q.time is for missing beat suggested by find_missing
    while (t <= to + sps) {
        if (countdown < 0 && sample_valid())
            qf();
        else if (countdown < 0) {
            // 1 second in samples.
            countdown = strtim("1");
            state = CLEANUP;
        }
        else if (countdown-- <= 0)
            break;


        q0 = q(t); q1 = q(t-1); q2 = q(t-2);

        // Found a peak
        if (q1 > pthr && q2 < q1 && q1 >= q0 && t > dt4) {
            // Add a peak if q1 is bigger than its neighbors and we are not
            // within dt4 of the signal start.
            addpeak(t-1, q1);
            last_peak = t-1;

            /*why does p start from cpeak->next? It gets around to the
            first peak but why not start from there?*/
            for (p = cpeak->next; p->time < t - rtmax; p = p->next) {

                /* If the peak is prominent and rrmin away from the previous
                qrs detection...

                (Limitation: Cannot detect qrs complexes at the start?
                unless searchback?)
                */
                if (p->time >= annot.time + rrmin && peaktype(p) == 1) {
                    // If the qrs-threshold is surpassed...
                    if (p->amp > qthr) {
                        rr = p->time - annot.time;

                        /* 
                            "it is suspected that a low-amplitude beat may have
                            been missed between two consecutive detected beats"
                        */
                        if (
                                rr > rrmean + 2 * rrdev && 
                                rr > 2 * (rrmean - rrdev) &&
                                (q = find_missing(r, p))) {
                            /* Revisited an old peak and judged as qrs*/
                            p = q;
                            rr = p->time - annot.time;
                            // Label it as a different type of beat found
                            annot.subtyp = 1;
                        }

                        // update rr mean and std
                        if ((rrd = rr - rrmean) < 0)
                            rrd = -rrd;
                        rrdev += (rrd - rrdev) >> 3;
                        
                        if (rrd > rrinc)
                            rrd = rrinc;
                        
                        if (rr > rrmean)
                            rrmean += rrd;
                        else
                            rrmean -= rrd;
                        
                        // Update qrs threshold
                        if (p->amp > qthr * 4)
                            qthr++;
                        else if (p->amp < qthr)
                            qthr--;
                        
                        if (qthr > pthr * 20)qthr = pthr * 20;

                        last_qrs = p->time;

                        // Write the annotation
                        if (state == RUNNING) {
                            int qsize;
                            annot.time = p->time - dt2;
                            qsize = p->amp * 10.0 / qthr;
                            if (qsize > 127)
                                qsize = 127;
                            annot.num = qsize;
                            putann(0, &annot);
                            annot.time += dt2;
                        }

                        /* look for this beat's T-wave (tw is a peak)*/
                        tw = NULL; rtdmin = rtmean;
                        for (q = p->next; q->time > annot.time; q = q->next) {
                            rt = q->time - annot.time - dt2;
                            if (rt < rtmin) continue;
                            if (rt > rtmax) break;
                            if ((rtd = rt - rtmean) < 0) rtd = -rtd;
                            if (rtd < rtdmin) {
                            rtdmin = rtd;
                            tw = q;
                            }
                        }
                        if (tw) {
                            static WFDB_Annotation tann;
                            tann.time = tw->time - dt2;
         
                            rt = tann.time - annot.time;
                            if ((rtmean += (rt - rtmean) >> 4) > rtmax)
                                rtmean = rtmax;
                            else if (rtmean < rtmin)
                                rtmean = rrmin;
                            tw->type = 2;   /* mark T-wave as secondary */
                        }

                        // r is set to the peak if its amplitude passes qthr
                        r = p; q = NULL; qamp = 0; annot.subtyp = 0;

                    }
                    // peak amplitude not above qthr.
                    // reduce q threshold if it's been longer than rrmax
                    else if (t - last_qrs > rrmax && qthr > qthmin){
                        qthr -= (qthr >> 4);
                    }
                }
            }
        }
        // sample is not a peak
        // if it's been too long since we found a peak and the peak threshold
        // is not at its min, reduce it by 4bit shift (like subtracting 1/16)
        else if (t - last_peak > rrmax && pthr > pthmin)
            pthr -= (pthr >> 4);

    }

    if (state == LEARNING)
    return;
    
    /* Mark the last beat or two. */
    for (p = cpeak->next; p->time < (p->next)->time; p = p->next) {
    //  if (to > 0 && p->time > to - sps)
    //    break;
    if (p->time >= annot.time + rrmin && p->time < tf && peaktype(p) == 1) {
        annot.anntyp = NORMAL;
        annot.chan = sig;
        annot.time = p->time;p
        putann(0, &annot);
    }
    }
}

