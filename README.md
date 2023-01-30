# cold_pool_variogram

## Description
Code to calculate empirical variogram and "gradiogram" (mean gradients) for temperature observations of FESSTVaL 2021 station network. Can also be applied to any regionalized variable at any given locations.

## Usage
```python
import pandas as pd
import empirical_variogram as eva

# Read station coordinates and create pairs
coordinates = eva.read_coordinates('network_coordinates.txt')
pairs       = eva.pair_stations(coordinates)

# Set up synthetic temperature data
temp_data = pd.DataFrame(columns=coordinates.index)
for t in range(0,5): temp_data.loc[t] = coordinates['X'] * (t+1) * (0.1/1000) + 20

# Calculate variogram and gradiogram (mean gradients) for single time step
temp_variogram  = eva.calc_variogram(temp_data.loc[0],pairs)
temp_gradiogram = eva.calc_variogram(temp_data.loc[0],pairs,func=eva.gradiogram)

# Calculate variogram for time period
temp_variogram_period = eva.calc_variogram_period(temp_data.loc[0:3],pairs)

```

## Contact
Bastian Kirsch (bastian.kirsch@uni-hamburg.de)<br>
Meteorologisches Institut, Universität Hamburg, Germany

30 January 2023
