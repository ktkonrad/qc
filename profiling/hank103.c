/* hank103.f -- translated by f2c (version 20090411).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static doublecomplex c_b2 = {2.,0.};
static integer c__9 = 9;
static integer c__2 = 2;
static integer c__18 = 18;

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */


/*        this is the end of the debugging code and the beginning of the */
/*        hankel function code proper. */


/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */





/* Subroutine */ int hanks103_(doublecomplex *z__, doublecomplex *hanks, 
	integer *n, integer *ifexpon)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, i1;
    static doublecomplex cd, cdd;
    extern /* Subroutine */ int hank103_(doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *);


/*       This subroutine evaluates the first n+1 Hankel functions of the */
/*       argument z. The user also has the option of evaluating the */
/*       functions H_m(z) scaled by the (complex) coefficient e^{-i \cdot z}. */
/*       This option is provided via the parameter ifexpon (see below) */


/*                      input parameters: */

/*  z - the complex number for which the hankel functions */
/*        H_0, H_1 are to be evaluated */
/*  n - the highest order of any Hankel function to be evaluated */
/*  ifexpon - the integer parameter telling the subroutine whether */
/*        to calculate the actual values of the hankel functions, */
/*        or the values of Hankel functions scaled by e^{-i \cdot z}. */
/*        Permitted values: 0 and 1. */
/*    ifexpon = 1 will cause the subroutine to evaluate the Hankel functions */
/*        honestly */
/*    ifexpon = 0 will cause the subroutine to scale the Hankel functions */
/*        by e^{-i \cdot z}. */

/*                      output parameters: */

/*  hanks - the first n+1 Hankel functions of the (complex) argument z. */
/*        Please note that hanks(1) is the Hankel function of order 0, */
/*        hanks(2) is the Hankel function of order 1, ..., hanks(n+1) */
/*        is the Hankel function of order n */

/*       . . . evaluate the functions h0,h1 */

    /* Parameter adjustments */
    --hanks;

    /* Function Body */
    hank103_(z__, &hanks[1], &hanks[2], ifexpon);


/*       conduct recursion */

    z_div(&z__1, &c_b2, z__);
    cd.r = z__1.r, cd.i = z__1.i;
    cdd.r = cd.r, cdd.i = cd.i;
    i__1 = *n;
    for (i1 = 2; i1 <= i__1; ++i1) {

	i__ = i1 - 1;

/* ccc        hanks(i1+1)=(2*i)/z*hanks(i1)-hanks(i1-1) */
	i__2 = i1 + 1;
	i__3 = i1;
	z__2.r = cdd.r * hanks[i__3].r - cdd.i * hanks[i__3].i, z__2.i = 
		cdd.r * hanks[i__3].i + cdd.i * hanks[i__3].r;
	i__4 = i1 - 1;
	z__1.r = z__2.r - hanks[i__4].r, z__1.i = z__2.i - hanks[i__4].i;
	hanks[i__2].r = z__1.r, hanks[i__2].i = z__1.i;

	z__1.r = cdd.r + cd.r, z__1.i = cdd.i + cd.i;
	cdd.r = z__1.r, cdd.i = z__1.i;
/* L1200: */
    }

    return 0;
} /* hanks103_ */






/* Subroutine */ int hank103_(doublecomplex *z__, doublecomplex *h0, 
	doublecomplex *h1, integer *ifexpon)
{
    /* Initialized data */

    static doublecomplex ima = {0.,1.};
    static doublereal pi = 3.1415926535897932;

    /* System generated locals */
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;
    static doublecomplex equiv_0[1];

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *), z_exp(doublecomplex *, 
	    doublecomplex *), z_div(doublecomplex *, doublecomplex *, 
	    doublecomplex *), z_log(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublecomplex y0, y1, z2, cd, zr, zu, fj0, fj1, h0r, h1r, h0u, h1u;
#define rea ((doublereal *)equiv_0)
#define com (equiv_0)
    static integer ier;
    static doublecomplex ser2, ser3;
    static doublereal half, subt;
    static doublecomplex cclog;
    extern /* Subroutine */ int hank103r_(doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *), hank103u_(
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *);


/*        this subroutine evaluates the hankel functions H_0^1, H_1^1 */
/*        for an arbitrary user-specified complex number z. The user */
/*        also has the option of evaluating the functions h0, h1 */
/*        scaled by the (complex) coefficient e^{-i \cdot z}. This */
/*        subroutine is a modification of the subroutine hank102 */
/*        (see), different from the latter by having the parameter */
/*        ifexpon. Please note that the subroutine hank102 is in */
/*        turn a slightly accelerated version of the old hank101 */
/*        (see). The principal claim to fame of all three is that */
/*        they are valid on the whole  complex plane, and are */
/*        reasonably accurate (14-digit relative accuracy) and */
/*        reasonably fast. Also, please note that all three have not */
/*        been carefully tested in the third quadrant (both x and y */
/*        negative); some sort of numerical trouble is possible */
/*        (though has not been observed) for LARGE z in the third */
/*        quadrant. */

/*                      input parameters: */

/*  z - the complex number for which the hankel functions */
/*        H_0, H_1 are to be evaluated */
/*  ifexpon - the integer parameter telling the subroutine whether */
/*        to calculate the actual values of the hankel functions, */
/*        or the values of Hankel functions scaled by e^{-i \cdot z}. */
/*        Permitted values: 0 and 1. */
/*    ifexpon = 1 will cause the subroutine to evaluate the Hankel functions */
/*        honestly */
/*    ifexpon = 0 will cause the subroutine to scale the Hankel functions */
/*        by e^{-i \cdot z}. */

/*                      output parameters: */

/*  h0, h1 - the said Hankel functions */


/*        . . . if z in the upper half-plane - act accordingly */

    com->r = z__->r, com->i = z__->i;
    if (rea[1] < 0.) {
	goto L1400;
    }
    hank103u_(z__, &ier, h0, h1, ifexpon);
    return 0;
L1400:

/*       if z is in the right lower quadrant - act accordingly */

    if (rea[0] < 0.) {
	goto L2000;
    }
    hank103r_(z__, &ier, h0, h1, ifexpon);
    return 0;
L2000:

/*       z is in the left lower quadrant. compute */
/*       h0, h1 at the points zu, zr obtained from z by reflection */
/*       in the x and y axis, respectively */

    d_cnjg(&z__1, z__);
    zu.r = z__1.r, zu.i = z__1.i;
    z__1.r = -zu.r, z__1.i = -zu.i;
    zr.r = z__1.r, zr.i = z__1.i;

    hank103u_(&zu, &ier, &h0u, &h1u, ifexpon);
    hank103r_(&zr, &ier, &h0r, &h1r, ifexpon);
    if (*ifexpon == 1) {
	goto L3000;
    }
    com->r = zu.r, com->i = zu.i;
    subt = abs(rea[1]);
    z__3.r = ima.r * zu.r - ima.i * zu.i, z__3.i = ima.r * zu.i + ima.i * 
	    zu.r;
    z__2.r = z__3.r - subt, z__2.i = z__3.i;
    z_exp(&z__1, &z__2);
    cd.r = z__1.r, cd.i = z__1.i;
    z__1.r = h0u.r * cd.r - h0u.i * cd.i, z__1.i = h0u.r * cd.i + h0u.i * 
	    cd.r;
    h0u.r = z__1.r, h0u.i = z__1.i;
    z__1.r = h1u.r * cd.r - h1u.i * cd.i, z__1.i = h1u.r * cd.i + h1u.i * 
	    cd.r;
    h1u.r = z__1.r, h1u.i = z__1.i;
    z__3.r = ima.r * zr.r - ima.i * zr.i, z__3.i = ima.r * zr.i + ima.i * 
	    zr.r;
    z__2.r = z__3.r - subt, z__2.i = z__3.i;
    z_exp(&z__1, &z__2);
    cd.r = z__1.r, cd.i = z__1.i;
    z__1.r = h0r.r * cd.r - h0r.i * cd.i, z__1.i = h0r.r * cd.i + h0r.i * 
	    cd.r;
    h0r.r = z__1.r, h0r.i = z__1.i;
    z__1.r = h1r.r * cd.r - h1r.i * cd.i, z__1.i = h1r.r * cd.i + h1r.i * 
	    cd.r;
    h1r.r = z__1.r, h1r.i = z__1.i;
L3000:

/*       compute the functions j0, j1, y0, y1 */
/*       at the point zr */

    half = 1.;
    half /= 2;
    z__3.r = h0u.r + h0r.r, z__3.i = h0u.i + h0r.i;
    z__2.r = half * z__3.r, z__2.i = half * z__3.i;
    z_div(&z__1, &z__2, &ima);
    y0.r = z__1.r, y0.i = z__1.i;
    z__3.r = h0u.r - h0r.r, z__3.i = h0u.i - h0r.i;
    z__2.r = -z__3.r, z__2.i = -z__3.i;
    z__1.r = half * z__2.r, z__1.i = half * z__2.i;
    fj0.r = z__1.r, fj0.i = z__1.i;

    z__4.r = h1u.r - h1r.r, z__4.i = h1u.i - h1r.i;
    z__3.r = -z__4.r, z__3.i = -z__4.i;
    z__2.r = half * z__3.r, z__2.i = half * z__3.i;
    z_div(&z__1, &z__2, &ima);
    y1.r = z__1.r, y1.i = z__1.i;
    z__2.r = h1u.r + h1r.r, z__2.i = h1u.i + h1r.i;
    z__1.r = half * z__2.r, z__1.i = half * z__2.i;
    fj1.r = z__1.r, fj1.i = z__1.i;

/*        finally, compute h0, h1 */

/*       . . . calculate ser2, ser3 */

    d_cnjg(&z__2, z__);
    z__1.r = -z__2.r, z__1.i = -z__2.i;
    z2.r = z__1.r, z2.i = z__1.i;
    z_log(&z__1, &z2);
    cclog.r = z__1.r, cclog.i = z__1.i;
    d__1 = 2.;
    z__4.r = d__1 * fj0.r, z__4.i = d__1 * fj0.i;
    z__3.r = z__4.r / pi, z__3.i = z__4.i / pi;
    z__2.r = z__3.r * cclog.r - z__3.i * cclog.i, z__2.i = z__3.r * cclog.i + 
	    z__3.i * cclog.r;
    z__1.r = y0.r - z__2.r, z__1.i = y0.i - z__2.i;
    ser2.r = z__1.r, ser2.i = z__1.i;
    d__1 = 2.;
    z__4.r = d__1 * fj1.r, z__4.i = d__1 * fj1.i;
    z__3.r = z__4.r / pi, z__3.i = z__4.i / pi;
    z__2.r = z__3.r * cclog.r - z__3.i * cclog.i, z__2.i = z__3.r * cclog.i + 
	    z__3.i * cclog.r;
    z__1.r = y1.r - z__2.r, z__1.i = y1.i - z__2.i;
    ser3.r = z__1.r, ser3.i = z__1.i;

/*       reflect all of these in the imaginary axis */

    d_cnjg(&z__1, &fj0);
    fj0.r = z__1.r, fj0.i = z__1.i;
    d_cnjg(&z__2, &fj1);
    z__1.r = -z__2.r, z__1.i = -z__2.i;
    fj1.r = z__1.r, fj1.i = z__1.i;

    d_cnjg(&z__1, &ser2);
    ser2.r = z__1.r, ser2.i = z__1.i;
    d_cnjg(&z__2, &ser3);
    z__1.r = -z__2.r, z__1.i = -z__2.i;
    ser3.r = z__1.r, ser3.i = z__1.i;

/*       reconstitute y0, y1 */

    z_log(&z__1, z__);
    cclog.r = z__1.r, cclog.i = z__1.i;
    d__1 = 2.;
    z__4.r = d__1 * fj0.r, z__4.i = d__1 * fj0.i;
    z__3.r = z__4.r / pi, z__3.i = z__4.i / pi;
    z__2.r = z__3.r * cclog.r - z__3.i * cclog.i, z__2.i = z__3.r * cclog.i + 
	    z__3.i * cclog.r;
    z__1.r = ser2.r + z__2.r, z__1.i = ser2.i + z__2.i;
    y0.r = z__1.r, y0.i = z__1.i;
    d__1 = 2.;
    z__4.r = d__1 * fj1.r, z__4.i = d__1 * fj1.i;
    z__3.r = z__4.r / pi, z__3.i = z__4.i / pi;
    z__2.r = z__3.r * cclog.r - z__3.i * cclog.i, z__2.i = z__3.r * cclog.i + 
	    z__3.i * cclog.r;
    z__1.r = ser3.r + z__2.r, z__1.i = ser3.i + z__2.i;
    y1.r = z__1.r, y1.i = z__1.i;

    z__2.r = ima.r * y0.r - ima.i * y0.i, z__2.i = ima.r * y0.i + ima.i * 
	    y0.r;
    z__1.r = fj0.r + z__2.r, z__1.i = fj0.i + z__2.i;
    h0->r = z__1.r, h0->i = z__1.i;
    z__2.r = ima.r * y1.r - ima.i * y1.i, z__2.i = ima.r * y1.i + ima.i * 
	    y1.r;
    z__1.r = fj1.r + z__2.r, z__1.i = fj1.i + z__2.i;
    h1->r = z__1.r, h1->i = z__1.i;
    if (*ifexpon == 1) {
	return 0;
    }
    z__4.r = -ima.r, z__4.i = -ima.i;
    z__3.r = z__4.r * z__->r - z__4.i * z__->i, z__3.i = z__4.r * z__->i + 
	    z__4.i * z__->r;
    z__2.r = z__3.r + subt, z__2.i = z__3.i;
    z_exp(&z__1, &z__2);
    cd.r = z__1.r, cd.i = z__1.i;
    z__1.r = h0->r * cd.r - h0->i * cd.i, z__1.i = h0->r * cd.i + h0->i * 
	    cd.r;
    h0->r = z__1.r, h0->i = z__1.i;
    z__1.r = h1->r * cd.r - h1->i * cd.i, z__1.i = h1->r * cd.i + h1->i * 
	    cd.r;
    h1->r = z__1.r, h1->i = z__1.i;
    return 0;
} /* hank103_ */

#undef com
#undef rea







/* Subroutine */ int hank103u_(doublecomplex *z__, integer *ier, 
	doublecomplex *h0, doublecomplex *h1, integer *ifexpon)
{
    /* Initialized data */

    static doublecomplex ima = {0.,1.};
    static struct {
	doublereal e_1[70];
	} equiv_1 = { -6.619836118357782e-13, -6.619836118612709e-13, 
		-7.3075142647542e-22, 3.928160926261892e-11, 
		5.712712520172854e-10, -5.712712519967086e-10, 
		-1.083820384008718e-8, -1.894529309455499e-19, 
		7.528123700585197e-8, 7.528123700841491e-8, 
		1.356544045548053e-17, -8.147940452202855e-7, 
		-3.568198575016769e-6, 3.568198574899888e-6, 
		2.592083111345422e-5, 4.2090748700194e-16, 
		-7.935843289157352e-5, -7.935843289415642e-5, 
		-6.848330800445365e-15, 4.136028298630129e-4, 
		9.210433149997867e-4, -9.210433149680665e-4, 
		-.003495306809056563, -6.469844672213905e-14, 
		.005573890502766937, .005573890503000873, 
		3.76734185797815e-13, -.01439178509436339, 
		-.01342403524448708, .01342403524340215, .008733016209933828, 
		1.400653553627576e-12, .02987361261932706, .02987361261607835,
		 -3.388096836339433e-12, -.1690673895793793, 
		.2838366762606121, -.2838366762542546, .7045107746587499, 
		-5.363893133864181e-12, -.7788044738211666, -.778804473813036,
		 5.524779104964783e-12, 1.146003459721775, .6930697486173089, 
		-.6930697486240221, -.7218270272305891, 3.633022466839301e-12,
		 .3280924142354455, .3280924142319602, -1.472323059106612e-12,
		 -.2608421334424268, -.09031397649230536, .09031397649339185, 
		.05401342784296321, -3.464095071668884e-13, 
		-.01377057052946721, -.01377057052927901, 
		4.273263742980154e-14, .005877224130705015, 
		.001022508471962664, -.001022508471978459, 
		-2.789107903871137e-4, 2.283984571396129e-15, 
		2.799719727019427e-5, 2.7997197269709e-5, 
		-3.371218242141487e-17, -3.682310515545645e-6, 
		-1.191412910090512e-7, 1.191412910113518e-7 };

    static struct {
	doublereal e_1[70];
	} equiv_4 = { 4.428361927253983e-13, -4.428361927153559e-13, 
		-2.575693161635231e-11, -2.878656317479645e-22, 
		3.658696304107867e-10, 3.658696304188925e-10, 
		7.463138750413651e-20, -6.748894854135266e-9, 
		-4.530098210372099e-8, 4.530098210271137e-8, 
		4.698787882823243e-7, 5.343848349451927e-18, 
		-1.948662942158171e-6, -1.948662942204214e-6, 
		-1.658085463182409e-16, 1.31690610049657e-5, 
		3.645368564036497e-5, -3.645368563934748e-5, 
		-1.63345854781839e-4, -2.697770638600506e-15, 
		2.81678497655166e-4, 2.816784976676616e-4, 
		2.54867335118006e-14, -6.106478245116582e-4, 
		2.054057459296899e-4, -2.054057460218446e-4, 
		-.00625496236729126, 1.484073406594994e-13, 
		.01952900562500057, .01952900562457318, 
		-5.517611343746895e-13, -.08528074392467523, 
		-.1495138141086974, .1495138141099772, .4394907314508377, 
		-1.334677126491326e-12, -1.113740586940341, 
		-1.113740586937837, 2.113005088866033e-12, 1.170212831401968, 
		1.262152242318805, -1.262152242322008, -1.557810619605511, 
		2.176383208521897e-12, .8560741701626648, .8560741701600203, 
		-1.431161194996653e-12, -.8386735092525187, -.365181917659929,
		 .3651819176613019, .2811692367666517, -5.799941348040361e-13,
		 -.0949463018293728, -.0949463018289448, 
		1.364615527772751e-13, .05564896498129176, .01395239688792536,
		 -.0139523968879995, -.005871314703753967, 
		1.683372473682212e-14, .001009157100083457, 
		.001009157100077235, -8.997331160162008e-16, 
		-2.723724213360371e-4, -2.708696587599713e-5, 
		2.70869658761883e-5, 3.533092798326666e-6, 
		-1.328028586935163e-17, -1.134616446885126e-7, 
		-1.134616446876064e-7 };

    static struct {
	doublereal e_1[62];
	} equiv_6 = { .5641895835516786, -.564189583551601, 
		-3.902447089770041e-10, -3.334441074447365e-12, 
		-.07052368835911731, -.07052368821797083, 1.95729931508537e-9,
		 -3.126801711815631e-7, -.03967331737107949, 
		.03967327747706934, 6.902866639752817e-5, 
		3.178420816292497e-7, .0408045716606128, .04080045784614144, 
		-2.218731025620065e-5, .006518438331871517, 
		.09798339748600499, -.09778028374972253, -.3151825524811773, 
		-7.995603166188139e-4, 1.111323666639636, 1.11679117899433, 
		.01635711249533488, -8.527067497983841, -25.95553689471247, 
		25.86942834408207, 134.5583522428299, .2002017907999571, 
		-308.6364384881525, -309.4609382885628, -1.505974589617013, 
		1250.150715797207, 2205.210257679573, -2200.328091885836, 
		-6724.941072552172, -7.018887749450317, 8873.498980910335, 
		8891.369384353965, 20.08805099643591, -20306.81426035686, 
		-20100.17782384992, 20060.46282661137, 34279.41581102808, 
		34.32892927181724, -25114.17407338804, -25165.67363193558, 
		-33.18253740485142, 31439.40826027085, 16584.66564673543, 
		-16548.43151976437, -14463.4504132651, -16.45433213663233, 
		5094.709396573681, 5106.816671258367, 3.470692471612145, 
		-2797.902324245621, -561.5581955514127, 560.1021281020627, 
		146.3856702925587, .1990076422327786, -9.334741618922085, 
		-9.361368967669095 };

    static struct {
	doublereal e_1[62];
	} equiv_8 = { -.5641895835446003, -.5641895835437973, 
		3.473016376419171e-11, -3.710264617214559e-10, 
		.2115710836381847, -.2115710851180242, 3.132928887334847e-7, 
		2.064187785625558e-8, -.06611954881267806, -.0661199717690031,
		 -3.38600489318156e-6, 7.146557892862998e-5, 
		-.05728505088320786, .05732906930408979, -.006884187195973806,
		 -2.383737409286457e-4, .1170452203794729, .1192356405185651, 
		.008652871239920498, -.3366165876561572, -1.203989383538728, 
		1.144625888281483, 9.153684260534125, .1781426600949249, 
		-27.40411284066946, -28.34461441294877, -2.19261107160634, 
		144.5470231392735, 336.1116314072906, -327.0584743216529, 
		-1339.254798224146, -16.57618537130453, 2327.097844591252, 
		2380.960024514808, 77.60611776965994, -7162.513471480693, 
		-9520.608696419367, 9322.604506839242, 21440.33447577134, 
		223.0232555182369, -20875.84364240919, -21317.62020653283, 
		-382.5699231499171, 35829.76792594737, 26426.32405857713, 
		-25851.37938787267, -32514.46505037506, -371.0875194432116, 
		16838.05377643986, 17243.93921722052, 184.6128226280221, 
		-14797.35877145448, -5258.288893282565, 5122.237462705988, 
		2831.540486197358, 39.05972651440027, -556.2781548969544, 
		-572.6891190727206, -2.246192560136119, 146.5347141877978, 
		9.456733342595993, -9.155767836700837 };


    /* System generated locals */
    doublecomplex z__1, z__2, z__3;
    static doublecomplex equiv_2[1];

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *), z_sqrt(doublecomplex *, 
	    doublecomplex *), z_div(doublecomplex *, doublecomplex *, 
	    doublecomplex *), z_exp(doublecomplex *, doublecomplex *), pow_zi(
	    doublecomplex *, doublecomplex *, integer *);

    /* Local variables */
    static doublereal d__;
    static integer m;
    static doublecomplex cd;
#define c0p1 ((doublereal *)&equiv_1)
#define c1p1 ((doublereal *)&equiv_4)
#define c0p2 ((doublereal *)&equiv_6)
#define rea ((doublereal *)equiv_2)
#define c1p2 ((doublereal *)&equiv_8)
#define com (equiv_2)
#define c0p1b ((doublereal *)&equiv_1 + 34)
#define c1p1b ((doublereal *)&equiv_4 + 34)
#define c0p2b ((doublereal *)&equiv_6 + 34)
#define c1p2b ((doublereal *)&equiv_8 + 34)
#define buf01 ((doublereal *)&equiv_1 + 33)
#define buf11 ((doublereal *)&equiv_4 + 33)
#define buf02 ((doublereal *)&equiv_6 + 33)
#define buf12 ((doublereal *)&equiv_8 + 33)
    static doublecomplex ccex;
    static doublereal done;
    static doublecomplex zzz9;
    extern /* Subroutine */ int hank103a_(doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *), hank103l_(doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *), hank103p_(
	    doublereal *, integer *, doublecomplex *, doublecomplex *);
    static doublereal thresh1, thresh2, thresh3;


/*        this subroutine evaluates the hankel functions H_0^1, H_1^1 */
/*        for a user-specified complex number z in the upper half-plane. */
/*        it is reasonably accurate (14-digit relative accuracy) */
/*        and reasonably fast. */


/*                      input parameters: */

/*  z - the complex number for which the hankel functions */
/*        H_0, H_1 are to be evaluated */

/*                      output parameters: */

/*  ier - error return code. */
/*         ier=0 means successful conclusion */
/*         ier=4 means that z is not in the upper half-plane */
/*  h0, h1 - the said Hankel functions */






/*        if the user-specified z is in the lower half-plane */
/*        - bomb out */

    *ier = 0;
    com->r = z__->r, com->i = z__->i;
    if (rea[1] >= 0.) {
	goto L1200;
    }
    *ier = 4;
    return 0;
L1200:

    done = 1.;
    thresh1 = 1.;
    thresh2 = 13.690000000000001f;
    thresh3 = 400.;

/*       check if if the user-specified z is in one of the */
/*       intermediate regimes */

    d_cnjg(&z__2, z__);
    z__1.r = z__->r * z__2.r - z__->i * z__2.i, z__1.i = z__->r * z__2.i + 
	    z__->i * z__2.r;
    d__ = z__1.r;
    if (d__ < thresh1 || d__ > thresh3) {
	goto L3000;
    }

/*        the user-specified z is in one of the intermediate regimes. */
/*        act accordingly */


    if (d__ > thresh2) {
	goto L2000;
    }

/*       z is in the first intermediate regime: its absolute value is */
/*       between 1 and 3.7. act accordingly */

/*       . . . evaluate the expansion */

    z__2.r = done, z__2.i = 0.;
    z_sqrt(&z__3, z__);
    z_div(&z__1, &z__2, &z__3);
    cd.r = z__1.r, cd.i = z__1.i;

    ccex.r = cd.r, ccex.i = cd.i;
    if (*ifexpon == 1) {
	z__3.r = ima.r * z__->r - ima.i * z__->i, z__3.i = ima.r * z__->i + 
		ima.i * z__->r;
	z_exp(&z__2, &z__3);
	z__1.r = ccex.r * z__2.r - ccex.i * z__2.i, z__1.i = ccex.r * z__2.i 
		+ ccex.i * z__2.r;
	ccex.r = z__1.r, ccex.i = z__1.i;
    }

    pow_zi(&z__1, z__, &c__9);
    zzz9.r = z__1.r, zzz9.i = z__1.i;
    m = 35;
    hank103p_(c0p1, &m, &cd, h0);
    z__2.r = h0->r * ccex.r - h0->i * ccex.i, z__2.i = h0->r * ccex.i + h0->i 
	    * ccex.r;
    z__1.r = z__2.r * zzz9.r - z__2.i * zzz9.i, z__1.i = z__2.r * zzz9.i + 
	    z__2.i * zzz9.r;
    h0->r = z__1.r, h0->i = z__1.i;

    hank103p_(c1p1, &m, &cd, h1);
    z__2.r = h1->r * ccex.r - h1->i * ccex.i, z__2.i = h1->r * ccex.i + h1->i 
	    * ccex.r;
    z__1.r = z__2.r * zzz9.r - z__2.i * zzz9.i, z__1.i = z__2.r * zzz9.i + 
	    z__2.i * zzz9.r;
    h1->r = z__1.r, h1->i = z__1.i;
    return 0;
L2000:

/*       z is in the second intermediate regime: its absolute value is */
/*       between 3.7 and 20. act accordingly. */

    z__2.r = done, z__2.i = 0.;
    z_sqrt(&z__3, z__);
    z_div(&z__1, &z__2, &z__3);
    cd.r = z__1.r, cd.i = z__1.i;

    ccex.r = cd.r, ccex.i = cd.i;
    if (*ifexpon == 1) {
	z__3.r = ima.r * z__->r - ima.i * z__->i, z__3.i = ima.r * z__->i + 
		ima.i * z__->r;
	z_exp(&z__2, &z__3);
	z__1.r = ccex.r * z__2.r - ccex.i * z__2.i, z__1.i = ccex.r * z__2.i 
		+ ccex.i * z__2.r;
	ccex.r = z__1.r, ccex.i = z__1.i;
    }
    m = 31;
    hank103p_(c0p2, &m, &cd, h0);
    z__1.r = h0->r * ccex.r - h0->i * ccex.i, z__1.i = h0->r * ccex.i + h0->i 
	    * ccex.r;
    h0->r = z__1.r, h0->i = z__1.i;

    m = 31;
    hank103p_(c1p2, &m, &cd, h1);
    z__1.r = h1->r * ccex.r - h1->i * ccex.i, z__1.i = h1->r * ccex.i + h1->i 
	    * ccex.r;
    h1->r = z__1.r, h1->i = z__1.i;
    return 0;
L3000:

/*        z is either in the local regime or the asymptotic one. */
/*        if it is in the local regime - act accordingly. */

    if (d__ > 50.) {
	goto L4000;
    }
    hank103l_(z__, h0, h1, ifexpon);
    return 0;

/*        z is in the asymptotic regime. act accordingly. */

L4000:
    hank103a_(z__, h0, h1, ifexpon);
    return 0;
} /* hank103u_ */

#undef buf12
#undef buf02
#undef buf11
#undef buf01
#undef c1p2b
#undef c0p2b
#undef c1p1b
#undef c0p1b
#undef com
#undef c1p2
#undef rea
#undef c0p2
#undef c1p1
#undef c0p1






/* Subroutine */ int hank103p_(doublecomplex *p, integer *m, doublecomplex *
	z__, doublecomplex *f)
{
    /* System generated locals */
    integer i__1;
    doublecomplex z__1, z__2;

    /* Local variables */
    static integer i__;


/*       evaluate a polynomial at a point */

    /* Parameter adjustments */
    --p;

    /* Function Body */
    i__1 = *m;
    f->r = p[i__1].r, f->i = p[i__1].i;
    for (i__ = *m - 1; i__ >= 1; --i__) {
	z__2.r = f->r * z__->r - f->i * z__->i, z__2.i = f->r * z__->i + f->i 
		* z__->r;
	i__1 = i__;
	z__1.r = z__2.r + p[i__1].r, z__1.i = z__2.i + p[i__1].i;
	f->r = z__1.r, f->i = z__1.i;
/* L1200: */
    }
    return 0;
} /* hank103p_ */






/* Subroutine */ int hank103a_(doublecomplex *z__, doublecomplex *h0, 
	doublecomplex *h1, integer *ifexpon)
{
    /* Initialized data */

    static doublecomplex ima = {0.,1.};
    static doublereal pi = 3.1415926535897932;
    static doublereal done = 1.;
    static doublecomplex cdumb = {.70710678118654757,-.70710678118654746};
    static doublereal p[18] = { 1.,-.0703125,.112152099609375,
	    -.5725014209747314,6.074042001273483,-110.0171402692467,
	    3038.090510922384,-118838.4262567833,6252951.493434797,
	    -425939216.5047669,36468400807.06556,-3833534661393.944,
	    485401468685290.1,-72868573493776570.,1.279721941975975e19,
	    -2.599382102726235e21,6.046711487532401e23,-1.597065525294211e26 }
	    ;
    static doublereal q[18] = { -.125,.0732421875,-.2271080017089844,
	    1.727727502584457,-24.38052969955606,551.3358961220206,
	    -18257.75547429317,832859.3040162893,-50069589.53198893,
	    3836255180.230434,-364901081884.9834,42189715702840.96,
	    -5827244631566907.,9.47628809926011e17,-1.792162323051699e20,
	    3.900121292034e22,-9.677028801069847e24,2.715581773544907e27 };
    static doublereal p1[18] = { 1.,.1171875,-.144195556640625,
	    .6765925884246826,-6.883914268109947,121.5978918765359,
	    -3302.272294480852,127641.2726461746,-6656367.718817687,
	    450278600.3050393,-38338575207.42789,4011838599133.198,
	    -506056850331472.6,75726164611179570.,-1.326257285320556e19,
	    2.687496750276277e21,-6.2386705823747e23,1.644739123064188e26 };
    static doublereal q1[18] = { .375,-.1025390625,.2775764465332031,
	    -1.993531733751297,27.24882731126854,-603.8440767050702,
	    19718.37591223663,-890297.8767070679,53104110.10968522,
	    -4043620325.107754,382701134659.8606,-44064814178522.79,
	    6065091351222699.,-9.83388387659068e17,1.855045211579829e20,
	    -4.027994121281017e22,9.974783533410457e24,-2.794294288720121e27 }
	    ;

    /* System generated locals */
    integer i__1;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4, z__5;
    static doublecomplex equiv_0[1];

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), pow_zi(
	    doublecomplex *, doublecomplex *, integer *), z_exp(doublecomplex 
	    *, doublecomplex *), z_sqrt(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, m;
    static doublecomplex pp, qq, pp1, qq1, cdd;
#define rea ((doublereal *)equiv_0)
#define com (equiv_0)
    static doublecomplex zinv, zinv22, cccexp;





/*        evaluate the asymptotic expansion for h0,h1 at */
/*        the user-supplied point z, provided it is not */
/*        in the fourth quadrant */

    m = 10;
    z__2.r = done, z__2.i = 0.;
    z_div(&z__1, &z__2, z__);
    zinv.r = z__1.r, zinv.i = z__1.i;

    i__1 = m - 1;
    pp.r = p[i__1], pp.i = 0.;
    i__1 = m - 1;
    pp1.r = p1[i__1], pp1.i = 0.;
    pow_zi(&z__1, &zinv, &c__2);
    zinv22.r = z__1.r, zinv22.i = z__1.i;

    i__1 = m - 1;
    qq.r = q[i__1], qq.i = 0.;
    i__1 = m - 1;
    qq1.r = q1[i__1], qq1.i = 0.;

    for (i__ = m - 1; i__ >= 1; --i__) {
	z__2.r = pp.r * zinv22.r - pp.i * zinv22.i, z__2.i = pp.r * zinv22.i 
		+ pp.i * zinv22.r;
	i__1 = i__ - 1;
	z__1.r = z__2.r + p[i__1], z__1.i = z__2.i;
	pp.r = z__1.r, pp.i = z__1.i;
	z__2.r = pp1.r * zinv22.r - pp1.i * zinv22.i, z__2.i = pp1.r * 
		zinv22.i + pp1.i * zinv22.r;
	i__1 = i__ - 1;
	z__1.r = z__2.r + p1[i__1], z__1.i = z__2.i;
	pp1.r = z__1.r, pp1.i = z__1.i;
	z__2.r = qq.r * zinv22.r - qq.i * zinv22.i, z__2.i = qq.r * zinv22.i 
		+ qq.i * zinv22.r;
	i__1 = i__ - 1;
	z__1.r = z__2.r + q[i__1], z__1.i = z__2.i;
	qq.r = z__1.r, qq.i = z__1.i;
	z__2.r = qq1.r * zinv22.r - qq1.i * zinv22.i, z__2.i = qq1.r * 
		zinv22.i + qq1.i * zinv22.r;
	i__1 = i__ - 1;
	z__1.r = z__2.r + q1[i__1], z__1.i = z__2.i;
	qq1.r = z__1.r, qq1.i = z__1.i;
/* L1600: */
    }

    z__1.r = qq.r * zinv.r - qq.i * zinv.i, z__1.i = qq.r * zinv.i + qq.i * 
	    zinv.r;
    qq.r = z__1.r, qq.i = z__1.i;
    z__1.r = qq1.r * zinv.r - qq1.i * zinv.i, z__1.i = qq1.r * zinv.i + qq1.i 
	    * zinv.r;
    qq1.r = z__1.r, qq1.i = z__1.i;

    cccexp.r = 1., cccexp.i = 0.;
    if (*ifexpon == 1) {
	z__2.r = ima.r * z__->r - ima.i * z__->i, z__2.i = ima.r * z__->i + 
		ima.i * z__->r;
	z_exp(&z__1, &z__2);
	cccexp.r = z__1.r, cccexp.i = z__1.i;
    }

    d__1 = 2 / pi;
    z__2.r = d__1 * zinv.r, z__2.i = d__1 * zinv.i;
    z_sqrt(&z__1, &z__2);
    cdd.r = z__1.r, cdd.i = z__1.i;

    z__2.r = ima.r * qq.r - ima.i * qq.i, z__2.i = ima.r * qq.i + ima.i * 
	    qq.r;
    z__1.r = pp.r + z__2.r, z__1.i = pp.i + z__2.i;
    h0->r = z__1.r, h0->i = z__1.i;
    z__3.r = cdd.r * cdumb.r - cdd.i * cdumb.i, z__3.i = cdd.r * cdumb.i + 
	    cdd.i * cdumb.r;
    z__2.r = z__3.r * cccexp.r - z__3.i * cccexp.i, z__2.i = z__3.r * 
	    cccexp.i + z__3.i * cccexp.r;
    z__1.r = z__2.r * h0->r - z__2.i * h0->i, z__1.i = z__2.r * h0->i + 
	    z__2.i * h0->r;
    h0->r = z__1.r, h0->i = z__1.i;

    z__2.r = ima.r * qq1.r - ima.i * qq1.i, z__2.i = ima.r * qq1.i + ima.i * 
	    qq1.r;
    z__1.r = pp1.r + z__2.r, z__1.i = pp1.i + z__2.i;
    h1->r = z__1.r, h1->i = z__1.i;
    z__5.r = -cdd.r, z__5.i = -cdd.i;
    z__4.r = z__5.r * cccexp.r - z__5.i * cccexp.i, z__4.i = z__5.r * 
	    cccexp.i + z__5.i * cccexp.r;
    z__3.r = z__4.r * cdumb.r - z__4.i * cdumb.i, z__3.i = z__4.r * cdumb.i + 
	    z__4.i * cdumb.r;
    z__2.r = z__3.r * h1->r - z__3.i * h1->i, z__2.i = z__3.r * h1->i + 
	    z__3.i * h1->r;
    z__1.r = z__2.r * ima.r - z__2.i * ima.i, z__1.i = z__2.r * ima.i + 
	    z__2.i * ima.r;
    h1->r = z__1.r, h1->i = z__1.i;

    return 0;
} /* hank103a_ */

#undef com
#undef rea







/* Subroutine */ int hank103l_(doublecomplex *z__, doublecomplex *h0, 
	doublecomplex *h1, integer *ifexpon)
{
    /* Initialized data */

    static doublereal gamma = .5772156649015328606;
    static doublecomplex ima = {0.,1.};
    static doublereal pi = 3.1415926535897932;
    static doublereal two = 2.;
    static doublereal cj0[16] = { 1.,-.25,.015625,-4.340277777777778e-4,
	    6.781684027777778e-6,-6.781684027777778e-8,4.709502797067901e-10,
	    -2.402807549524439e-12,9.385966990329841e-15,
	    -2.896903392077112e-17,7.242258480192779e-20,
	    -1.496334396734045e-22,2.597802772107717e-25,
	    -3.842903509035085e-28,4.901662639075363e-31,
	    -5.446291821194848e-34 };
    static doublereal cj1[16] = { -.5,.0625,-.002604166666666667,
	    5.425347222222222e-5,-6.781684027777778e-7,5.651403356481481e-9,
	    -3.363930569334215e-11,1.501754718452775e-13,
	    -5.214426105738801e-16,1.448451696038556e-18,
	    -3.291935672814899e-21,6.234726653058522e-24,
	    -9.991549123491221e-27,1.372465538941102e-29,
	    -1.633887546358454e-32,1.70196619412339e-35 };
    static doublereal ser2[16] = { .25,-.0234375,7.957175925925926e-4,
	    -1.41285083912037e-5,1.548484519675926e-7,-1.153828185281636e-9,
	    6.230136717695511e-12,-2.550971742728932e-14,
	    8.195247730999099e-17,-2.121234517551702e-19,
	    4.518746345057852e-22,-8.06152930228997e-25,1.222094716680443e-27,
	    -1.593806157473552e-30,1.807204342667468e-33,
	    -1.798089518115172e-36 };
    static doublereal ser2der[16] = { .5,-.09375,.004774305555555556,
	    -1.130280671296296e-4,1.548484519675926e-6,-1.384593822337963e-8,
	    8.722191404773715e-11,-4.081554788366291e-13,
	    1.475144591579838e-15,-4.242469035103405e-18,
	    9.941241959127275e-21,-1.934767032549593e-23,
	    3.177446263369152e-26,-4.462657240925946e-29,
	    5.421613028002404e-32,-5.75388645796855e-35 };

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    void pow_zi(doublecomplex *, doublecomplex *, integer *), z_log(
	    doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *), z_exp(doublecomplex *, 
	    doublecomplex *);

    /* Local variables */
    static integer i__, m;
    static doublecomplex y0, y1, z2, cd, fj0, fj1, cdddlog;



/*        this subroutine evaluates the hankel functions H_0^1, H_1^1 */
/*        for a user-specified complex number z in the local regime, */
/*        i. e. for cdabs(z) < 1 in the upper half-plane, */
/*        and for cdabs(z) < 4 in the lower half-plane, */
/*        it is reasonably accurate (14-digit relative accuracy) and */
/*        reasonably fast. */

/*                      input parameters: */

/*  z - the complex number for which the hankel functions */
/*        H_0, H_1 are to be evaluated */

/*                      output parameters: */

/*  h0, h1 - the said Hankel functions */


/*        evaluate j0, j1 */

    m = 16;
    fj0.r = 0., fj0.i = 0.;
    fj1.r = 0., fj1.i = 0.;
    y0.r = 0., y0.i = 0.;
    y1.r = 0., y1.i = 0.;
    pow_zi(&z__1, z__, &c__2);
    z2.r = z__1.r, z2.i = z__1.i;
    cd.r = 1., cd.i = 0.;

    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	z__2.r = cj0[i__2] * cd.r, z__2.i = cj0[i__2] * cd.i;
	z__1.r = fj0.r + z__2.r, z__1.i = fj0.i + z__2.i;
	fj0.r = z__1.r, fj0.i = z__1.i;
	i__2 = i__ - 1;
	z__2.r = cj1[i__2] * cd.r, z__2.i = cj1[i__2] * cd.i;
	z__1.r = fj1.r + z__2.r, z__1.i = fj1.i + z__2.i;
	fj1.r = z__1.r, fj1.i = z__1.i;
	i__2 = i__ - 1;
	z__2.r = ser2der[i__2] * cd.r, z__2.i = ser2der[i__2] * cd.i;
	z__1.r = y1.r + z__2.r, z__1.i = y1.i + z__2.i;
	y1.r = z__1.r, y1.i = z__1.i;
	z__1.r = cd.r * z2.r - cd.i * z2.i, z__1.i = cd.r * z2.i + cd.i * 
		z2.r;
	cd.r = z__1.r, cd.i = z__1.i;
	i__2 = i__ - 1;
	z__2.r = ser2[i__2] * cd.r, z__2.i = ser2[i__2] * cd.i;
	z__1.r = y0.r + z__2.r, z__1.i = y0.i + z__2.i;
	y0.r = z__1.r, y0.i = z__1.i;
/* L1800: */
    }
    z__2.r = -fj1.r, z__2.i = -fj1.i;
    z__1.r = z__2.r * z__->r - z__2.i * z__->i, z__1.i = z__2.r * z__->i + 
	    z__2.i * z__->r;
    fj1.r = z__1.r, fj1.i = z__1.i;

    z__3.r = z__->r / two, z__3.i = z__->i / two;
    z_log(&z__2, &z__3);
    z__1.r = z__2.r + gamma, z__1.i = z__2.i;
    cdddlog.r = z__1.r, cdddlog.i = z__1.i;
    z__2.r = cdddlog.r * fj0.r - cdddlog.i * fj0.i, z__2.i = cdddlog.r * 
	    fj0.i + cdddlog.i * fj0.r;
    z__1.r = z__2.r + y0.r, z__1.i = z__2.i + y0.i;
    y0.r = z__1.r, y0.i = z__1.i;
    d__1 = two / pi;
    z__1.r = d__1 * y0.r, z__1.i = d__1 * y0.i;
    y0.r = z__1.r, y0.i = z__1.i;

    z__1.r = y1.r * z__->r - y1.i * z__->i, z__1.i = y1.r * z__->i + y1.i * 
	    z__->r;
    y1.r = z__1.r, y1.i = z__1.i;

    z__4.r = -cdddlog.r, z__4.i = -cdddlog.i;
    z__3.r = z__4.r * fj1.r - z__4.i * fj1.i, z__3.i = z__4.r * fj1.i + 
	    z__4.i * fj1.r;
    z_div(&z__5, &fj0, z__);
    z__2.r = z__3.r + z__5.r, z__2.i = z__3.i + z__5.i;
    z__1.r = z__2.r + y1.r, z__1.i = z__2.i + y1.i;
    y1.r = z__1.r, y1.i = z__1.i;
    z__3.r = -y1.r, z__3.i = -y1.i;
    z__2.r = two * z__3.r, z__2.i = two * z__3.i;
    z__1.r = z__2.r / pi, z__1.i = z__2.i / pi;
    y1.r = z__1.r, y1.i = z__1.i;

    z__2.r = ima.r * y0.r - ima.i * y0.i, z__2.i = ima.r * y0.i + ima.i * 
	    y0.r;
    z__1.r = fj0.r + z__2.r, z__1.i = fj0.i + z__2.i;
    h0->r = z__1.r, h0->i = z__1.i;
    z__2.r = ima.r * y1.r - ima.i * y1.i, z__2.i = ima.r * y1.i + ima.i * 
	    y1.r;
    z__1.r = fj1.r + z__2.r, z__1.i = fj1.i + z__2.i;
    h1->r = z__1.r, h1->i = z__1.i;

    if (*ifexpon == 1) {
	return 0;
    }

    z__3.r = -ima.r, z__3.i = -ima.i;
    z__2.r = z__3.r * z__->r - z__3.i * z__->i, z__2.i = z__3.r * z__->i + 
	    z__3.i * z__->r;
    z_exp(&z__1, &z__2);
    cd.r = z__1.r, cd.i = z__1.i;
    z__1.r = h0->r * cd.r - h0->i * cd.i, z__1.i = h0->r * cd.i + h0->i * 
	    cd.r;
    h0->r = z__1.r, h0->i = z__1.i;
    z__1.r = h1->r * cd.r - h1->i * cd.i, z__1.i = h1->r * cd.i + h1->i * 
	    cd.r;
    h1->r = z__1.r, h1->i = z__1.i;

    return 0;
} /* hank103l_ */






/* Subroutine */ int hank103r_(doublecomplex *z__, integer *ier, 
	doublecomplex *h0, doublecomplex *h1, integer *ifexpon)
{
    /* Initialized data */

    static doublecomplex ima = {0.,1.};
    static struct {
	doublereal e_1[70];
	} equiv_1 = { -4.268441995428495e-24, 4.374027848105921e-24, 
		9.876152216238049e-24, -1.065264808278614e-21, 
		6.240598085551175e-20, 6.65852998549011e-20, 
		-5.107210870050163e-18, -2.931746613593983e-19, 
		1.611018217758854e-16, -1.359809022054077e-16, 
		-7.718746693707326e-16, 6.759496139812828e-15, 
		-1.067620915195442e-13, -1.434699000145826e-13, 
		3.868453040754264e-12, 7.06185339258518e-13, 
		-6.220133527871203e-11, 3.957226744337817e-11, 
		3.080863675628417e-10, -1.1546184312819e-9, 
		7.793319486868695e-9, 1.502570745460228e-8, 
		-1.97809085263843e-7, -7.39669187349903e-8, 
		2.175857247417038e-6, -8.473534855334919e-7, 
		-1.05338132760972e-5, 2.042555121261223e-5, 
		-4.812568848956982e-5, -1.961519090873697e-4, 
		.001291714391689374, 9.23442238495005e-4, -.01113890671502769,
		 9.053687375483149e-4, .05030666896877862, 
		-.04923119348218356, .5202355973926321, -.1705244841954454, 
		-1.134990486611273, -1.747542851820576, 8.308174484970718, 
		2.952358687641577, -32.86074510100263, 11.26542966971545, 
		65.76015458463394, -100.6116996293757, 32.16834899377392, 
		361.4005342307463, -665.3878500833375, -688.3582242804924, 
		2193.362007156572, 242.3724600546293, -3665.925878308203, 
		2474.933189642588, 1987.663383445796, -7382.586600895061, 
		4991.253411017503, 10085.05017740918, -12852.84928905621, 
		-5153.67482166847, 13016.56757246985, -4821.250366504323, 
		-4982.112643422311, 9694.070195648748, -1685.723189234701, 
		-6065.143678129265, 2029.510635584355, 1244.402339119502, 
		-433.6682903961364, 89.23209875101459 };

    static struct {
	doublereal e_1[70];
	} equiv_4 = { -4.019450270734195e-24, -4.819240943285824e-24, 
		1.087220822839791e-21, 1.219058342725899e-22, 
		-7.458149572694168e-20, 5.677825613414602e-20, 
		8.351815799518541e-19, -5.188585543982425e-18, 
		1.221075065755962e-16, 1.789261470637227e-16, 
		-6.829972121890858e-15, -1.497462301804588e-15, 
		1.579028042950957e-13, -9.4149603037588e-14, 
		-1.127570848999746e-12, 3.883137940932639e-12, 
		-3.397569083776586e-11, -6.779059427459179e-11, 
		1.149529442506273e-9, 4.363087909873751e-10, 
		-1.620182360840298e-8, 6.404695607668289e-9, 
		9.651461037419628e-8, -1.948572160668177e-7, 
		6.397881896749446e-7, 2.318661930507743e-6, 
		-1.983192412396578e-5, -1.294811208715315e-5, 
		2.062663873080766e-4, -2.867633324735777e-5, 
		-.001084309075952914, .001227880935969686, 
		2.538406015667726e-4, -.01153316815955356, .04520140008266983,
		 .05693944718258218, -.9640790976658534, -.6517135574036008, 
		2.051491829570049, -1.124151010077572, -3.977380460328048, 
		8.200665483661009, -7.950131652215817, -35.03037697046647, 
		96.07320812492044, 78.9407968985807, -374.9002890488298, 
		-8.153831134140778, 782.4282518763973, -603.5276543352174, 
		-500.4685759675768, 2219.009060854551, -2111.301101664672, 
		-4035.632271617418, 7319.737262526823, 2878.734389521922, 
		-10874.04934318719, 3945.740567322783, 6727.823761148537, 
		-12535.55346597302, 3440.468371829973, 13832.40926370073, 
		-9324.927373036743, -6181.580304530313, 6376.198146666679, 
		-1033.615527971958, -1497.604891055181, 1929.025541588262, 
		-42.19760183545219, -452.1162915353207 };

    static struct {
	doublereal e_1[54];
	} equiv_6 = { .5641895835569398, -.5641895835321127, 
		-.07052370223565544, -.07052369923405479, -.03966909368581382,
		 .03966934297088857, .04130698137268744, .04136196771522681, 
		.06240742346896508, -.06553556513852438, -.03258849904760676, 
		-.07998036854222177, -3.98800631195527, 1.327373751674479, 
		61.21789346915312, -92.51865216627577, 424.7064992018806, 
		2692.55333348915, -43746.91601489926, -36252.48208112831, 
		1010975.818048476, -28593.60062580096, -11389702.41206912, 
		10510979.79526042, 22840388.99211195, -203801251.5235694, 
		1325194353.842857, 1937443530.361381, -22459990186.52171, 
		-5998903865.344352, 179323705487.6609, -86251598823.06147, 
		-588776304273.5203, 1345331284205.28, -2743432269370.813, 
		-8894942160272.255, 42764631137945.64, 26650198866477.81, 
		-228072742395549.8, 36869087905539.73, 563984631816861.5, 
		-684152905161570.3, 99014267999660.38, 2798406605978152., 
		-4910062244008171., -5126937967581805., 13872929519367560., 
		1043295727224325., -15652041206872650., 12152628069735770., 
		3133802397107054., -18013945508070780., 4427598668012807., 
		6923499968336864. };

    static struct {
	doublereal e_1[62];
	} equiv_8 = { -.564189583543198, -.5641895835508094, 
		.2115710934750869, -.2115710923186134, -.06611607335011594, 
		-.06611615414079688, -.05783289433408652, .05785737744023628, 
		.08018419623822896, .08189816020440689, .1821045296781145, 
		-.217973897300874, .5544705668143094, 2.22446631644444, 
		-85.63271248520645, -43.94325758429441, 2720.62754707134, 
		-670.5390850875292, -39362.2196060077, 57917.30432605451, 
		-197678.7738827811, -1502498.631245144, 21553178.23990686, 
		18709537.96705298, -470399571.1098311, 3716595.90645319, 
		5080557859.012385, -4534199223.888966, -10644382116.47413, 
		86122438937.45942, -546601768778.5078, -807095038664.0701, 
		9337074941225.827, 2458379240643.264, -75486921712445.79, 
		37510931699543.36, 246067743135003.9, -599191937288191.1, 
		1425679408434606., 4132221939781502., -22475064694689690., 
		-12697710781650260., 129733629274902600., -28026269097913080.,
		 -346713722281301700., 477395521558219200., 
		-234716577658020600., -2.233638097535785e18, 
		5.382350866778548e18, 4.820328886922998e18, 
		-1.928978948099345e19, 157549874775090700., 
		3.049162180215152e19, -2.837046201123502e19, 
		-5.429391644354291e18, 6.974653380104308e19, 
		-5.322120857794536e19, -6.739879079691706e19, 
		6.780343087166473e19, 1.053455984204666e19, 
		-2.218784058435737e19, 1.505391868530062e19 };


    /* System generated locals */
    doublecomplex z__1, z__2, z__3;
    static doublecomplex equiv_2[1];

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *), z_exp(doublecomplex *, 
	    doublecomplex *), z_sqrt(doublecomplex *, doublecomplex *), z_div(
	    doublecomplex *, doublecomplex *, doublecomplex *), pow_zi(
	    doublecomplex *, doublecomplex *, integer *);

    /* Local variables */
    static doublereal d__;
    static integer m;
    static doublecomplex cd, cdd;
#define c0p1 ((doublereal *)&equiv_1)
#define c1p1 ((doublereal *)&equiv_4)
#define c0p2 ((doublereal *)&equiv_6)
#define rea ((doublereal *)equiv_2)
#define c1p2 ((doublereal *)&equiv_8)
#define com (equiv_2)
    static doublecomplex zz18;
#define c0p1b ((doublereal *)&equiv_1 + 34)
#define c1p1b ((doublereal *)&equiv_4 + 34)
#define c0p2b ((doublereal *)&equiv_6 + 34)
#define c1p2b ((doublereal *)&equiv_8 + 34)
#define buf01 ((doublereal *)&equiv_1 + 33)
#define buf11 ((doublereal *)&equiv_4 + 33)
#define buf02 ((doublereal *)&equiv_6 + 33)
#define buf12 ((doublereal *)&equiv_8 + 33)
    static doublereal done;
    static doublecomplex cccexp;
    extern /* Subroutine */ int hank103a_(doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *), hank103l_(doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *), hank103p_(
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *);
    static doublereal thresh1, thresh2, thresh3;


/*        this subroutine evaluates the hankel functions H_0^1, H_1^1 */
/*        for a user-specified complex number z in the right lower */
/*        quadrant. it is reasonably accurate (14-digit relative */
/*        accuracy) and reasonably fast. */


/*                      input parameters: */

/*  z - the complex number for which the hankel functions */
/*        H_0, H_1 are to be evaluated */

/*                      output parameters: */

/*  ier - error return code. */
/*         ier=0 means successful conclusion */
/*         ier=4 means that z is not in the right lower quadrant */
/*  h0, h1 - the said Hankel functions */









/*        if z is not in the right lower quadrant - bomb out */

    *ier = 0;
    com->r = z__->r, com->i = z__->i;
    if (rea[0] >= 0. && rea[1] <= 0.) {
	goto L1400;
    }
    *ier = 4;
    return 0;
L1400:

    done = 1.;
    thresh1 = 16.;
    thresh2 = 64.;
    thresh3 = 400.;

/*       check if if the user-specified z is in one of the */
/*       intermediate regimes */

    d_cnjg(&z__2, z__);
    z__1.r = z__->r * z__2.r - z__->i * z__2.i, z__1.i = z__->r * z__2.i + 
	    z__->i * z__2.r;
    d__ = z__1.r;
    if (d__ < thresh1 || d__ > thresh3) {
	goto L3000;
    }

/*        if the user-specified z is in the first intermediate regime */
/*        (i.e. if its absolute value is between 4 and 8), act accordingly */

    if (d__ > thresh2) {
	goto L2000;
    }

    cccexp.r = 1., cccexp.i = 0.;
    if (*ifexpon == 1) {
	z__2.r = ima.r * z__->r - ima.i * z__->i, z__2.i = ima.r * z__->i + 
		ima.i * z__->r;
	z_exp(&z__1, &z__2);
	cccexp.r = z__1.r, cccexp.i = z__1.i;
    }
    z__2.r = done, z__2.i = 0.;
    z_sqrt(&z__3, z__);
    z_div(&z__1, &z__2, &z__3);
    cdd.r = z__1.r, cdd.i = z__1.i;
    z__2.r = done, z__2.i = 0.;
    z_div(&z__1, &z__2, z__);
    cd.r = z__1.r, cd.i = z__1.i;
    pow_zi(&z__1, z__, &c__18);
    zz18.r = z__1.r, zz18.i = z__1.i;
    m = 35;
    hank103p_((doublecomplex*)c0p1, &m, &cd, h0);
    z__3.r = h0->r * cdd.r - h0->i * cdd.i, z__3.i = h0->r * cdd.i + h0->i * 
	    cdd.r;
    z__2.r = z__3.r * cccexp.r - z__3.i * cccexp.i, z__2.i = z__3.r * 
	    cccexp.i + z__3.i * cccexp.r;
    z__1.r = z__2.r * zz18.r - z__2.i * zz18.i, z__1.i = z__2.r * zz18.i + 
	    z__2.i * zz18.r;
    h0->r = z__1.r, h0->i = z__1.i;

    hank103p_((doublecomplex*)c1p1, &m, &cd, h1);
    z__3.r = h1->r * cdd.r - h1->i * cdd.i, z__3.i = h1->r * cdd.i + h1->i * 
	    cdd.r;
    z__2.r = z__3.r * cccexp.r - z__3.i * cccexp.i, z__2.i = z__3.r * 
	    cccexp.i + z__3.i * cccexp.r;
    z__1.r = z__2.r * zz18.r - z__2.i * zz18.i, z__1.i = z__2.r * zz18.i + 
	    z__2.i * zz18.r;
    h1->r = z__1.r, h1->i = z__1.i;
    return 0;
L2000:

/*       z is in the second intermediate regime (i.e. its */
/*       absolute value is between 8 and 20). act accordingly. */

    z__2.r = done, z__2.i = 0.;
    z_div(&z__1, &z__2, z__);
    cd.r = z__1.r, cd.i = z__1.i;
    z_sqrt(&z__1, &cd);
    cdd.r = z__1.r, cdd.i = z__1.i;
    cccexp.r = 1., cccexp.i = 0.;
    if (*ifexpon == 1) {
	z__2.r = ima.r * z__->r - ima.i * z__->i, z__2.i = ima.r * z__->i + 
		ima.i * z__->r;
	z_exp(&z__1, &z__2);
	cccexp.r = z__1.r, cccexp.i = z__1.i;
    }
    m = 27;

    hank103p_((doublecomplex*)c0p2, &m, &cd, h0);
    z__2.r = h0->r * cccexp.r - h0->i * cccexp.i, z__2.i = h0->r * cccexp.i + 
	    h0->i * cccexp.r;
    z__1.r = z__2.r * cdd.r - z__2.i * cdd.i, z__1.i = z__2.r * cdd.i + 
	    z__2.i * cdd.r;
    h0->r = z__1.r, h0->i = z__1.i;

    m = 31;
    hank103p_((doublecomplex*)c1p2, &m, &cd, h1);
    z__2.r = h1->r * cccexp.r - h1->i * cccexp.i, z__2.i = h1->r * cccexp.i + 
	    h1->i * cccexp.r;
    z__1.r = z__2.r * cdd.r - z__2.i * cdd.i, z__1.i = z__2.r * cdd.i + 
	    z__2.i * cdd.r;
    h1->r = z__1.r, h1->i = z__1.i;
    return 0;
L3000:


/*        z is either in the local regime or the asymptotic one. */
/*        if it is in the local regime - act accordingly. */

    if (d__ > 50.) {
	goto L4000;
    }
    hank103l_(z__, h0, h1, ifexpon);
    return 0;

/*        z is in the asymptotic regime. act accordingly. */

L4000:
    hank103a_(z__, h0, h1, ifexpon);
    return 0;
} /* hank103r_ */

#undef buf12
#undef buf02
#undef buf11
#undef buf01
#undef c1p2b
#undef c0p2b
#undef c1p1b
#undef c0p1b
#undef com
#undef c1p2
#undef rea
#undef c0p2
#undef c1p1
#undef c0p1


