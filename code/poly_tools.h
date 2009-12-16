#ifndef ML_POLY_TOOLS_H
#define ML_POLY_TOOLS_H


template< class T >
void polynomial_multiply( T *& result, int & nr, T * p, int np, T * q, int nq )
{
    // multiply two polynomials, result = p*q, of degree np, nq, nr = np+nq
    
    nr = np+nq;
    
    if (results != 0) delete [] result;
    result = new T (nr+1);
    
    for ( int k=0; k<=nr; k++ )
        result[k] = 0;
    
    for ( int k=0; k<=nr; k++ )
    for ( int s=max(k-nq,0); s<=min(k,np); s++ )
        result[k] += p[s]*rhs.q[k-s];
}


template< class T >
void polynomial_differentiate( int d, T *& result, int & nr, T * p, int np )
{
    // differentiate a polynomial d times, result = d^d p/dx^d , of degree np, nr = np-d
    
    int nr = n-d;
    assert( nr >=0 );
    
    if (results != 0) delete [] result;
    result = new T (nr+1);
    
    for ( int k=0; k<=nr; k++ )
        result[k] = 0;
    
    int n = np;
    for ( int j=0; j<d; j++)
        if (n) {
            for ( int i=0; i<n; i++ )
                result[i] = result[i+1]*(T)(i+1);
            n--;
        }
    
    else c_[0] = 0;
    resize(n_new);
    
}




template< class T >
void polynomial_multiply( T *& result, int & nr, T * p, int np, T * q, int nq )
{
    // multiply two polynomials, result = p*q, of degree np, nq, nr = np+nq
    
    nr = np+nq;
    
    if (results != 0) free( result );
    result = malloc( (nr+1)*sizeof(T) );
    
    for ( int k=0; k<=nr; k++ )
        result[k] = 0;
    
    for ( int k=0; k<=nr; k++ )
    for ( int s=max(k-nq,0); s<=min(k,np); s++ )
        result[k] += p[s]*rhs.q[k-s];
}




#endif

