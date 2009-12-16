
#include <mathlib/math/std_math.h>

#include <mathlib/math/hdaf/hdaf.h>
#include <mathlib/math/transforms/fft.h>
#include <mathlib/math/random/ml_random.h>
#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>


#define mp_0 ((mp_real)(0.0))
#define mp_1 ((mp_real)(1.0))
#define mp_2 ((mp_real)(2.0))
#define mp_pi (atan(mp_1)*4.0)
#define mp_2pi (mp_pi*2.0)
#define mp_4pi2 (mp_2pi*mp_2pi)
#define mp_sqrt2 (sqrt(mp_2))
#define mp_iu ( complex<mp_real > (mp_0, mp_1 ) )


/*


import math
import pylab as p
import sys
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *

p.plot ( read_file( "frequency" ),  read_file( "error" ) , "*r", read_file( "frequency" ), read_file( "expected_error" ) )


p.plot ( read_file( "frequency" ), read_file( "response" ), read_file( "frequency" ), read_file( "measured_response" ), "*r" )

p.show()



*/



void kernel_extension ( double *& ker_ext, int n, double * kernel, int n_l, int n_r )
{
    //  i*n2+j
    
    if (ker_ext == 0) ker_ext = ml_alloc<double> (n);
    
    for (int k=0; k<n; k++)
        ker_ext[k] = 0;
    
    for (int k=-n_l; k<0; k++)
    {
        int l = k;
        while ( l < 0 ) l += n;
        ker_ext[l ] += kernel[k];
    }
    
    for (int k=0; k<=n_r; k++)
        ker_ext[ k%n ] += kernel[k];
    
}


void hdaf_diff_kernel ( double *& kernel, int & w, double L, int n, int m, double k, int d )
{
    mp_real sigma = sqrt( mp_2*m+ mp_1)/k;
    
    ml_poly<mp_real > P;
    make_hdaf_ml_poly(P,m );
    differentiate_hdaf_poly( P, d );
    
    w = (int)ceil(hdaf_truncate_point (1E-16, 4E-16, m, dble(sigma), d )/(L/n));
    
    if (kernel != 0) ml_free( kernel );
    kernel = ml_alloc<double > (2*w+1);
    double * p = &(kernel[w]);
    
    mp_real h = ((mp_real)L)/n;
    mp_real s = h/(mp_sqrt2*sigma);
    mp_real ss = s*s;
    mp_real f = pow(mp_sqrt2*sigma,-(d+1))*(h/n);
    
    for (int k=-w; k<=w; k++ )
        p[k] = dble( exp(-ss*(k*k)) *P( dble(s*k) ) *f );
}

void hdaf_diff_kernel ( complex<double > *& fft_ker, double L, int n, int m, double k, int d )
{
    double * kernel = 0;
    int w = 0;
    
    hdaf_diff_kernel ( kernel, w, L, n, m, k, d );
    
    double * p = &(kernel[w]);
    double * ker_ext = 0;
    
    kernel_extension ( ker_ext, n, p, w, w );
    
    fft( ker_ext, fft_ker, n );
    
    ml_free (ker_ext);
    ml_free (kernel);
}

double l2_error( double *& exact, double *& approx, int n )
{
    double sum1 = 0.0;
    for (int k=0; k<n; k++)
        sum1 += (approx[k]-exact[k])*(approx[k]-exact[k]);
    
    double sum2 = 0.0;
    for (int k=0; k<n; k++)
        sum2 += exact[k]*exact[k];
    
    return sqrt(sum1/sum2);
}


double rms( double *& x, int n )
{
    double sum = 0.0;
    for (int k=0; k<n; k++)
        sum += x[k]*x[k];
    
    return sqrt(sum/n);
}

int main()
{
    std_setup();
    
    mp::mp_init(100);
    
    int n = pow(2,12);
    double a = 0.0;
    double b = ml_2pi;
    
    int d = 2;
    int m = 4;
    double gamma = 0.5;
    
    vector<double > frequency;
    vector<double > error;
    vector<double > measured_response;
    vector<double > response;
    vector<double > expected_error;
    
    
//    for (int k=1; k<=n/2; k += max(n/100,1) )
    
    for (int k=1; k<=200; k += 1 )
    {
        double K = (ml_2pi/(b-a))*k;
        
        double * signal = ml_alloc<double > (n);
        double * diff_signal = ml_alloc<double > (n);
        
        for (int i=0; i<n; i++)
        {
            double x = ((b-a)*i)/n+a;
            signal[i] = cos(x*K);
            diff_signal[i] = -cos(x*K)*K*K;
        }
        
        complex<double > * fft_ker=0;
        hdaf_diff_kernel ( fft_ker, b-a, n, m, (ml_2pi/(b-a))*(gamma*n/2), d );
        
        double * signal_filtered=0;
        complex<double > * workspace=0;
        
        fft( signal, workspace, n );
        
        for (int i=0; i<n/2+1; i++)
            workspace[i] *= fft_ker[i];
        
        ifft( workspace, signal_filtered, n );
        
        
        //---------------
        
        cout << K << "\t" << l2_error( diff_signal, signal_filtered, n ) << endl;
        
        hdaf_delta_hat delta( m, sqrt(2*m+1)/( (ml_2pi/(b-a))*(gamma*n/2) ) );
        
        frequency.push_back( K );
        measured_response.push_back( rms(signal_filtered,n)/rms(diff_signal,n ) );
        response.push_back ( delta(K) );
        
        expected_error.push_back ( log10( fabs(delta(K)-1.0)+1E-18 ) );
        error.push_back( log10(fabs(l2_error( diff_signal, signal_filtered, n )) +1E-18) );
        
        ml_free(fft_ker);
        ml_free(signal);
        
    }
    
    
    sprintf( fname, "/workspace/output/temp/frequency" );
    output ( frequency,  fname );
    sprintf( fname, "/workspace/output/temp/measured_response" );
    output ( measured_response,  fname );
    sprintf( fname, "/workspace/output/temp/error" );
    output ( error,  fname );
    sprintf( fname, "/workspace/output/temp/expected_error" );
    output ( expected_error,  fname );
    sprintf( fname, "/workspace/output/temp/response" );
    output ( response,  fname );
    
    std_exit();
}




