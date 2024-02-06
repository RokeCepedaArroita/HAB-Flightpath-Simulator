# Monte Carlo Simulation of Helium Balloon Flight

Monte Carlo simulation of the flight of a light payload suspended from a helium balloon created by the implementation of statistical variations in the standard Tawhiri predictor [[1]](http://tawhiri.cusf.co.uk/en/latest/), made for Manchester Students for the Development and Exploration of Space (MANSEDS) [[2]](http://manseds.co.uk/) for the launch of Mancunian Balloonian [[3]](http://manseds.co.uk/balloon_project....).

The 1σ, 2σ, and 3σ ellipses delineate the area within which the flight will be found at a particular time with a probability of 68%, 95%, and 99.7%, respectively. The traveled distance and the horizontal (\(v_h\)) and vertical (\(v_v\)) speeds correspond to the values experienced by the most likely flight.

The simulation consists of 10,000 statistical sets that arise from setting tolerances in the neck lift, launch time, launch location, drag coefficients, and payload and balloon masses. Variability in the burst altitude is modeled combining experimental data [[4]](http://www.southampton.ac.uk/~as7/pub...) and manufacturer data [[5]](https://ukhas.org.uk/guides:balloon_data). The descent is modeled by using an experimental descent speed distribution and neglecting parachute glide, which will be implemented in future versions. Weather model inaccuracies and uncertainties in the predictor itself are assumed to be much smaller than the intrinsic experimental errors.

This gives a very good tool to understand the uncertainties associated with a particular flight and can be used to assess correlations between flight paths and final uncertainties. It also enables a complete treatment and propagation of errors, which is a key step in determining the feasibility of a particular launch. The program generates the following launch card:

### Example Simulation Video
[![Monte Carlo Simulation Example](https://www.youtube.com/watch?v=9tVPAtmtc7w)](https://www.youtube.com/watch?v=9tVPAtmtc7w)

### Input
```plaintext
11-Feb-2017 11:30:00 ± 3 min
LAUNCH ALTITUDE = 85 ± 10 m
BALLOON MASS = 1274.0 ± 0.5 g
PAYLOAD MASS = 1111.4 ± 8.8 g
DESCENT RATE = 5.45 ± 1.40 m/s
```

### Output
```plaintext
DATASET = 06-Feb 6:00 AM (±0)
MEAN HELIUM VOLUME = 3.333 m^3
NECK LIFT = 2146.99 ± 0.50 g
ASCENT RATE = 5.00 ± 0.06 m/s
BURST ALTITUDE = 33567 ± 618 m
FLIGHT DURATION = 02:36 ± 10 min
AVE LANDING ACCURACY = 7.7 km
LAT LANDING ACCURACY = 7.2 km
LON LANDING ACCURACY = 2.8 km
1-SIGMA SEARCH AREA = 139.4 km^2
DISTANCE FROM LAUNCH = 98 km
TRAVELLED BY BALLOON = 142 km
ASCENT HEADING SHIFT = 370 deg
DESCENT HEADING SHIFT = 356 deg
```

## Links
1. Tawhiri predictor: [http://tawhiri.cusf.co.uk/en/latest/](http://tawhiri.cusf.co.uk/en/latest/)
2. MANSEDS: [http://manseds.co.uk/](http://manseds.co.uk/)
3. Mancunian Balloonian: [http://manseds.co.uk/balloon_project....](http://manseds.co.uk/balloon_project....)
4. Burst variability: [http://www.southampton.ac.uk/~as7/pub...](http://www.southampton.ac.uk/~as7/pub...)
5. Balloon data: [https://ukhas.org.uk/guides:balloon_data](https://ukhas.org.uk/guides:balloon_data)
```
