#ifndef _StringTok_H
#define _StringTok_H

template<class T>
class StringTok
{
public:
  StringTok( const T seq, 
    typename T::size_type pos = 0 )
    : seq_( seq ) , pos_( pos ) { }

  T operator()( const T& delim );
  
private:
  const T seq_;
  typename T::size_type pos_;
};

template<class T>
T StringTok<T>::operator()
  ( const T& delim )
{
  T token;

  if( pos_ != T::npos )
    {
      // start of found token

      typename T::size_type first = 
	seq_.find_first_not_of( delim.c_str(), pos_ );

      if( first != T::npos )
	{
	  // length of found token
	  typename T::size_type num =
	    seq_.find_first_of( delim.c_str(), first ) - first;

	  // do all the work off to the side
	  token = seq_.substr( first, num );    

	  // done; now commit using nonthrowing operations only
	  pos_ = first+num;                            
	  if( pos_ != T::npos ) ++pos_;                
	  if( pos_ >= seq_.size() ) pos_ = T::npos; 
	}
    }
  
  return token;
}

#endif
