# Thank you to Paul Schlyter and his resource (http://www.stjarnhimlen.se/comp/ppcomp.html) for making this possible

import math
import numpy
import geocoder
from datetime import datetime
from typing import Tuple

KM_IN_1_AU = 149597870.7
MILES_IN_1_KM = 0.62137119223733
EARTH_RADIUS_KM = 6371
SECONDS_IN_DAY = 86400
SPEED_OF_LIGHT_KMH = 1079252848.8

# Finds Moon's position and returns it as (distance_earth_radii, right_ascension, declination)
def moon_position() -> Tuple[float, float, float]:
    day = datetime_since_2000_start()

    longitude_ascending_node = (125.1228 - 0.0529538083 * day) % 360
    ecliptic_inclination = 5.1454
    argument_of_periapsis = (318.0634 + 0.1643573223 * day) % 360
    mean_distance_sun = 60.2666
    eccentricity = 0.054900
    mean_anomaly = (115.3654 + 13.0649929509 * day) % 360

    sun_argument_of_periapsis = 282.9404 + (4.70935*10**-5) * day
    sun_mean_anomaly = (356.0470 + 0.9856002585 * day) % 360
    sun_mean_longitude = (sun_argument_of_periapsis + sun_mean_anomaly) % 360
    mean_longitude = (longitude_ascending_node + argument_of_periapsis + mean_anomaly) % 360
    mean_elongation = (mean_longitude - sun_mean_longitude) % 360
    latitude_argument = (mean_longitude - longitude_ascending_node) % 360

    longitude_perturbations = (
        -1.274 * deg_sin(mean_anomaly - mean_elongation*2),
        0.658 * deg_sin(mean_elongation*2),
        -0.186 * deg_sin(sun_mean_anomaly),
        -0.059 * deg_sin(mean_anomaly*2 - mean_elongation*2),
        -0.057 * deg_sin(mean_anomaly - mean_elongation*2 + sun_mean_anomaly),
        0.053 * deg_sin(mean_anomaly + mean_elongation*2),
        0.046 * deg_sin(mean_elongation*2 - sun_mean_anomaly),
        0.041 * deg_sin(mean_anomaly - sun_mean_anomaly),
        -0.035 * deg_sin(mean_elongation),
        -0.031 * deg_sin(mean_anomaly + sun_mean_anomaly),
        -0.015 * deg_sin(latitude_argument*2 - mean_elongation*2),
        0.011 * deg_sin(mean_anomaly - mean_elongation*4)
    )

    latitude_perturbations = (
        -0.173 * deg_sin(latitude_argument - mean_elongation*2),
        -0.055 * deg_sin(mean_anomaly - latitude_argument - mean_elongation*2),
        -0.046 * deg_sin(mean_anomaly + latitude_argument - mean_elongation*2),
        0.033 * deg_sin(latitude_argument + mean_elongation*2),
        0.017 * deg_sin(mean_anomaly*2 + latitude_argument)
    )

    distance_perturbations = (
        -0.58 * deg_cos(mean_anomaly - mean_elongation*2),
        -0.46 * deg_cos(mean_elongation*2)
    )

    total_longitude_perturbations = sum(longitude_perturbations)
    total_latitude_perturbations = sum(latitude_perturbations)
    total_distance_perturbations = sum(distance_perturbations)

    eccentric_anomaly_first_approx = mean_anomaly + ((180/math.pi) * eccentricity * deg_sin(mean_anomaly)) * (1 + eccentricity * deg_cos(mean_anomaly))
    eccentric_anomaly = eccentric_anomaly_first_approx - (eccentric_anomaly_first_approx - ((180/math.pi) * eccentricity * deg_sin(eccentric_anomaly_first_approx)) - mean_anomaly) / (1 - eccentricity * deg_cos(eccentric_anomaly_first_approx))

    x_lunar_orbit_plane = mean_distance_sun * (deg_cos(eccentric_anomaly) - eccentricity)
    y_lunar_orbit_plane = mean_distance_sun * math.sqrt(1 - eccentricity**2) * deg_sin(eccentric_anomaly)

    distance_earth_radii = math.sqrt(x_lunar_orbit_plane**2 + y_lunar_orbit_plane**2)
    true_anomaly = numpy.arctan2(y_lunar_orbit_plane, x_lunar_orbit_plane) * (180/math.pi) % 360

    x_ecliptic = distance_earth_radii * (deg_cos(longitude_ascending_node) * deg_cos(true_anomaly+argument_of_periapsis) - deg_sin(longitude_ascending_node) * deg_sin(true_anomaly+argument_of_periapsis) * deg_cos(ecliptic_inclination))
    y_ecliptic = distance_earth_radii * (deg_sin(longitude_ascending_node) * deg_cos(true_anomaly+argument_of_periapsis) + deg_cos(longitude_ascending_node) * deg_sin(true_anomaly+argument_of_periapsis) * deg_cos(ecliptic_inclination))
    z_ecliptic = distance_earth_radii * deg_sin(true_anomaly+argument_of_periapsis) * deg_sin(ecliptic_inclination)

    longitude = math.atan2(y_ecliptic, x_ecliptic) * (180/math.pi) % 360 + total_longitude_perturbations
    latitude = math.atan2(z_ecliptic, math.sqrt(x_ecliptic**2 + y_ecliptic**2)) * (180/math.pi) + total_latitude_perturbations
    distance_earth_radii = math.sqrt(x_ecliptic**2 + y_ecliptic**2 + z_ecliptic**2) + total_distance_perturbations

    ecliptic_obliquity = 23.4393 - (3.563*10**-7) * day

    x_ecliptic = deg_cos(longitude) * deg_cos(latitude)
    y_ecliptic = deg_sin(longitude) * deg_cos(latitude)
    z_ecliptic = deg_sin(latitude)

    x_equatorial = x_ecliptic
    y_equatorial = y_ecliptic * deg_cos(ecliptic_obliquity) - z_ecliptic * deg_sin(ecliptic_obliquity)
    z_equatorial = y_ecliptic * deg_sin(ecliptic_obliquity) + z_ecliptic * deg_cos(ecliptic_obliquity)

    right_ascension = math.atan2(y_equatorial, x_equatorial) * (180/math.pi) % 360
    declination = math.atan2(z_equatorial, math.sqrt(x_equatorial**2 + y_equatorial**2)) * (180/math.pi)

    return distance_earth_radii, right_ascension, declination

# Finds Sun's position and returns it as (distance_AU, right_ascension, declination)
def sun_position() -> Tuple[float, float, float]:
    day = datetime_since_2000_start()

    argument_of_periapsis = 282.9404 + (4.70935*10**-5) * day
    eccentricity = 0.016709 - (1.151*10**-9) * day
    mean_anomaly = (356.0470 + 0.9856002585 * day) % 360

    ecliptic_obliquity = 23.4393 - (3.563*10**-7) * day

    eccentric_anomaly = mean_anomaly + (180/math.pi) * eccentricity * deg_sin(mean_anomaly) * (1 + eccentricity * deg_cos(mean_anomaly))

    x_ecliptic_plane = deg_cos(eccentric_anomaly) - eccentricity
    y_ecliptic_plane = deg_sin(eccentric_anomaly) * math.sqrt(1 - (eccentricity**2))

    distance_AU = math.sqrt(x_ecliptic_plane**2 + y_ecliptic_plane**2)
    true_anomaly = numpy.arctan2(y_ecliptic_plane, x_ecliptic_plane) * (180/math.pi)

    longitude = (true_anomaly + argument_of_periapsis) % 360

    rectangular_x = distance_AU * deg_cos(longitude)
    rectangular_y = distance_AU * deg_sin(longitude)

    x_equatorial = rectangular_x
    y_equatorial = rectangular_y * deg_cos(ecliptic_obliquity)
    z_equatorial = rectangular_y * deg_sin(ecliptic_obliquity)

    right_ascension = math.atan2(y_equatorial, x_equatorial) * (180/math.pi)
    declination =  math.atan2(z_equatorial, math.sqrt(x_equatorial**2 + y_equatorial**2)) * (180/math.pi)

    return distance_AU, right_ascension, declination

# Finds Sun's local azimuth and altitude and returns it as (azimuth, altitude)
def sun_azimuth_altitude() -> Tuple[float, float]:
    day = datetime_since_2000_start()

    user = geocoder.ip("me")
    longitude = user.lng
    latitude = user.lat

    hours_into_day = (day - int(day)) * 24

    argument_of_periapsis = 282.9404 + (4.70935*10**-5) * day
    mean_anomaly = (356.0470 + 0.9856002585 * day) % 360
    mean_longitude = (argument_of_periapsis + mean_anomaly) % 360

    greenwich_mean_sidereal_time_midnight = ((mean_longitude + 180) % 360) / 15
    sidereal_time = greenwich_mean_sidereal_time_midnight + hours_into_day + (longitude / 15)

    distance_AU, right_ascension, declination = sun_position()

    hour_angle = (sidereal_time - (right_ascension / 15)) * 15

    x = deg_cos(hour_angle) * deg_cos(declination)
    y = deg_sin(hour_angle) * deg_cos(declination)
    z = deg_sin(declination)

    x_horizontal = x * deg_sin(latitude) - z * deg_cos(latitude)
    y_horizontal = y
    z_horizontal = x * deg_cos(latitude) + z * deg_sin(latitude)

    azimuth = math.atan2(y_horizontal, x_horizontal) * (180/math.pi) + 180
    altitude = math.atan2(z_horizontal, math.sqrt(x_horizontal**2 + y_horizontal**2)) * (180/math.pi)
    return azimuth, altitude

# Calculates time since start of 2000 in days as a decimal number
def datetime_since_2000_start() -> float:
    start_of_2000 = datetime(2000, 1, 1)
    current_datetime = datetime.utcnow()

    datetime_passed = current_datetime - start_of_2000
    seconds_passed = datetime_passed.total_seconds()

    datetime_since_2000_start = (seconds_passed / SECONDS_IN_DAY) + 1
    return datetime_since_2000_start

# Determines time in minutes and seconds for light to reach Earth from Moon/Sun
def time_light_reach_earth(distance_from_earth: float) -> str:
    time_reach_earth = distance_from_earth / SPEED_OF_LIGHT_KMH
    minutes = time_reach_earth * 60
    seconds = (minutes - int(minutes)) * 60

    formatted_time = str(int(minutes)) + " minutes " + str(round(seconds, 1)) + " seconds"
    return formatted_time

# Converts angle in degrees to hours/degrees, minutes, and seconds
def format_angle(degrees: float, type: str) -> str:
    if type == "hours":
        angle = (degrees * 24) / 360
        unit_1 = "h "
        unit_2 = "m "
        unit_3 = "s"
    elif type == "degrees":
        angle = degrees
        unit_1 = u"\N{DEGREE SIGN} "
        unit_2 = "' "
        unit_3 = "\""
    first_category = int(angle)
    minutes = int((angle - first_category) * 60)
    seconds = round((angle - first_category - (minutes / 60)) * 3600)
    formatted_angle = str(first_category).zfill(2) + unit_1 + str(abs(minutes)).zfill(2) + unit_2 + str(abs(seconds)).zfill(2) + unit_3
    return formatted_angle

# Formats and prints data on Moon/Sun
def print_data(body_name: str):
    line = "_" * 55

    print(line)
    print()

    if body_name == "moon":
        distance_earth_radii, right_ascension, declination = moon_position()
        distance_in_km = distance_earth_radii * EARTH_RADIUS_KM
    elif body_name == "sun":
        distance_AU, right_ascension, declination = sun_position()
        distance_in_km = distance_AU * KM_IN_1_AU

    distance_in_miles = distance_in_km * MILES_IN_1_KM
    hms_str = format_angle(right_ascension, "hours")
    dms_str = format_angle(declination, "degrees")

    print("Distance: " + f"{round(distance_in_km):,}" + " km (" + f"{round(distance_in_miles):,}" + " mi)")
    print("Time for light to reach Earth: " + str(time_light_reach_earth(distance_in_km)))
    print("Right ascension: " + hms_str + " (" + str(round(right_ascension, 4)) + " degrees)")
    print("Declination: " + dms_str + " (" + str(round(declination, 4)) + " degrees)")

    if body_name == "sun":
        azimuth, altitude = sun_azimuth_altitude()
        print("Local azimuth: " + str(round(azimuth, 1)) + " degrees")
        print("Local altitude: " + str(round(altitude, 1)) + " degrees")

    print(line)
    print()

# Takes sine of number in degrees
def deg_sin(num_in_degrees: float) -> float:
    return math.sin(math.radians(num_in_degrees))

# Takes cosine of number in degrees
def deg_cos(num_in_degrees: float) -> float:
    return math.cos(math.radians(num_in_degrees))

# Thanks to the following website for the ASCII art text generation http://patorjk.com/software/taag/#p=display&f=Graffiti&t=Type%20Something%20
print(r"""
 _      _                          _                                   _           _   _____        _
| |    (_)               /\       | |                                 (_)         | | |  __ \      | |
| |     ___   _____     /  \   ___| |_ _ __ ___  _ __   ___  _ __ ___  _  ___ __ _| | | |  | | __ _| |_ __ _
| |    | \ \ / / _ \   / /\ \ / __| __| '__/ _ \| '_ \ / _ \| '_ ` _ \| |/ __/ _` | | | |  | |/ _` | __/ _` |
| |____| |\ V /  __/  / ____ \\__ \ |_| | | (_) | | | | (_) | | | | | | | (_| (_| | | | |__| | (_| | || (_| |
|______|_| \_/ \___| /_/    \_\___/\__|_|  \___/|_| |_|\___/|_| |_| |_|_|\___\__,_|_| |_____/ \__,_|\__\__,_|

""")

# Controls user interaction with command line
while True:
    print("Type one of the following:")
    print("    moon (shows current moon data)")
    print("    sun (shows current sun data)")

    user_input = input()

    if user_input == "moon":
        print_data("moon")
    elif user_input == "sun":
        print_data("sun")
    else:
        continue
