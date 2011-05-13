

/* Useful routines */

double lower_limit(double *a, int M)
{
        int i;
        double limit = 1e+308;
        
        for (i=1;i<=M;++i)
                if (a[i]<limit)
                        limit = a[i];

        return limit;
}

double upper_limit(double *a, int M)
{
        int i;
        double limit = -1e+308;
        
        for (i=1;i<=M;++i)
                if (a[i]>limit)
                        limit = a[i];

        return limit;
}

float flower_limit(float *a, int M)
{
        int i;
        float limit = 1e+38;
        
        for (i=1;i<=M;++i)
                if (a[i]<limit)
                        limit = a[i];

        return limit;
}

float fupper_limit(float *a, int M)
{
        int i;
        float limit = -1e+38;
        
        for (i=1;i<=M;++i)
                if (a[i]>limit)
                        limit = a[i];

        return limit;
}


void endian_array_4byte(char *array, int n)
/* flips 4 byte arrays big<->little endian, n values, beginning where
 * ptr points. Processes 4n bytes in total.
 */
{
	char *ptr,temp;
	int i;

	ptr = array;
	for (i=0; i<n; ++i,ptr+=4) {
		temp = ptr[0];
		ptr[0] = ptr[3];
		ptr[3] = temp;
		temp = ptr[1];
		ptr[1] = ptr[2];
		ptr[2] = temp;
	}
}

int endian_4byte(void *in)
/* returns flipped 4 byte value big<->little endian, without destroying orig.
 * Incoming arg is a pointer, returns value. 99/10/15
 */
{
	char *inc,outc[4];

	inc = (char *)in;
	outc[0] = inc[3];
	outc[1] = inc[2];
	outc[2] = inc[1];
	outc[3] = inc[0];

	return *((int *)outc);
}

double endian_8byte(double *in)
/* returns flipped 8 byte double big<->little endian, without destroying orig.
 * Incoming arg is a pointer, returns value. 99/10/15.
 */
{
	char *inc,outc[8];

	inc = (char *)in;
	outc[0] = inc[7];
	outc[1] = inc[6];
	outc[2] = inc[5];
	outc[3] = inc[4];
	outc[4] = inc[3];
	outc[5] = inc[2];
	outc[6] = inc[1];
	outc[7] = inc[0];

	return *((double *)outc);
}
