# nenergy, averages the energy of a time series in the new format
{s = s + $1}
{t = t + $2}
END {printf( "energy = %10.3E  time = %10.3E\n",s/NR,t/NR)}
