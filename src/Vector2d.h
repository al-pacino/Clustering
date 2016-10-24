#pragma once

///////////////////////////////////////////////////////////////////////////////

template<typename NUMERIC_TYPE>
struct CVector2d {
	typedef NUMERIC_TYPE NumericType;
	typedef NumericType DistanceType;

	NumericType X;
	NumericType Y;

	explicit CVector2d( NumericType x = 0, NumericType y = 0 ) :
		X( x ), Y( y )
	{
	}

	DistanceType Distance( const CVector2d& point ) const
	{
		const DistanceType dx = X - point.X;
		const DistanceType dy = Y - point.Y;
		return sqrt( dx * dx + dy * dy );
	}
};

///////////////////////////////////////////////////////////////////////////////
