############################################################################################################################
                                # 1. sun_moon_sgb
############################################################################################################################

# sun_moon.py MicroPython Port of lunarmath.c
# Calculate sun and moon rise and set times for any date and location

# Licensing and copyright: see README.md

# Source "Astronomy on the Personal Computer" by Montenbruck and Pfleger
# ISBN 978-3-540-67221-0

# Port from C++ to MicroPython performed by Peter Hinch 2023.
# Withcontributions from Raul Kompaß and Marcus Mendenhall: see
# https://github.com/orgs/micropython/discussions/13075
# Raul Kompaß perfomed major simplification of the maths for deriving rise and
# set_times with improvements in precision with 32-bit floats.

# Moon phase now in separate module

## Changes made by Simcha Gershon Bohrer marked with ##

import time
from math import sin, cos, sqrt, fabs, atan, radians, floor, pi,      atan2, degrees, asin ##


LAT = 53.29756504536339  # Local defaults
LONG = -2.102811634540558

# MicroPython wanton epochs:
# time.gmtime(0)[0] = 1970 or 2000 depending on platform.
# On CPython:
# (date(2000, 1, 1) - date(1970, 1, 1)).days * 24*60*60 = 946684800
# (date(2000, 1, 1) - date(1970, 1, 1)).days = 10957

# Return time now in days relative to platform epoch.
# System time is set to local time, and MP has no concept of this. Hence
# time.time() returns secs since epoch 00:00:00 local time. If lto is local time
# offset to UTC, provided -12 < lto < 12, the effect of rounding ensures the
# right number of days for platform epoch at UTC.
def now_days() -> int:
    secs_per_day = 86400  # 24 * 3600
    t = RiSet.mtime()  # Machine time as int. Can be overridden for test.
    t -= t % secs_per_day  # Previous Midnight
    return round(t / secs_per_day)  # Days since datum


def quad(ym, yz, yp):
    # See Astronomy on the PC P48-49, plus contribution from Marcus Mendenhall
    # finds the parabola throuh the three points (-1,ym), (0,yz), (1, yp)
    # and returns the values of x where the parabola crosses zero
    # (roots of the quadratic)
    # and the number of roots (0, 1 or 2) within the interval [-1, 1]
    nz = 0
    a = 0.5 * (ym + yp) - yz
    b = 0.5 * (yp - ym)
    c = yz

    if a == 0:  ## מניעת חלוקה באפס כי פעמים רבות הייתה שגיאה בקוד בקוטב הצפוני
        return 0, 0, 0, 0 ##
        
    xe = -b / (2 * a)
    ye = (a * xe + b) * xe + c
    dis = b * b - 4.0 * a * c  # discriminant of y=a*x^2 +bx +c
    if dis > 0:  # parabola has roots
        if b < 0:
            z2 = (-b + sqrt(dis)) / (2 * a)  # z2 is larger root in magnitude
        else:
            z2 = (-b - sqrt(dis)) / (2 * a)
        z1 = (c / a) / z2  # z1 is always closer to zero
        if fabs(z1) <= 1.0:
            nz += 1
        if fabs(z2) <= 1.0:
            nz += 1
        if z1 < -1.0:
            z1 = z2
        return nz, z1, z2, ye
    return 0, 0, 0, 0  # No roots


# **** GET MODIFIED JULIAN DATE FOR DAY RELATIVE TO TODAY ****

# Returns modified julian day number defined as mjd = jd - 2400000.5
# Deals only in integer MJD's: the JD of just after midnight will always end in 0.5
# hence the MJD of an integer day number will always be an integer

# Re platform comparisons get_mjd returns the same value regardless of
# the platform's epoch: integer days since 00:00 UTC on 17 November 1858.
def get_mjd(ndays: int = 0) -> int:
    secs_per_day = 86400  # 24 * 3600
    days_from_epoch = now_days() + ndays  # Days since platform epoch
    mjepoch = 40587  # Modified Julian date of C epoch (1 Jan 70)
    if time.gmtime(0)[0] == 2000:
        mjepoch += 10957
    return mjepoch + days_from_epoch  # Convert days from 1 Jan 70 to MJD


def frac(x):
    return x % 1


# Convert rise or set time to int. These can be None (no event).
def to_int(x):
    return None if x is None else round(x)


def minisun(t):
    # Output sin(dec), cos(dec), ra
    # returns the ra and dec of the Sun
    # in decimal hours, degs referred to the equinox of date and using
    # obliquity of the ecliptic at J2000.0 (small error for +- 100 yrs)
    # takes t centuries since J2000.0. Claimed good to 1 arcmin
    coseps = 0.9174805004
    sineps = 0.397780757938
    e = 0.0167  # אקסצנטריות מסלול הארץ סביב השמש ###########

    m = 2 * pi * frac(0.993133 + 99.997361 * t)
    dl = 6893.0 * sin(m) + 72.0 * sin(2 * m)
    l = 2 * pi * frac(0.7859453 + m / (2 * pi) + (6191.2 * t + dl) / 1296000)
    sl = sin(l)
    x = cos(l)
    y = coseps * sl
    z = sineps * sl
    # rho = sqrt(1 - z * z)
    # dec = (360.0 / 2 * pi) * atan(z / rho)
    # ra = ((48.0 / (2 * pi)) * atan(y / (x + rho))) % 24
    
    ############
    #distance_km = 149597870.7  # יחידה אסטרונומית אחת (AU) מרחק ממוצע בלבד
    # חישוב מרחק מהשמש לפי חוק קפלר
    distance_AU = (1 - e**2) / (1 + e * cos(m))
    distance_km = distance_AU * 149597870.7
    ############
    
    return x, y, z, distance_km ################## המרחק זו תוספת שלי


def minimoon(t):
    # takes t and returns the geocentric ra and dec
    # claimed good to 5' (angle) in ra and 1' in dec
    # tallies with another approximate method and with ICE for a couple of dates
    arc = 206264.8062
    coseps = 0.9174805004
    sineps = 0.397780757938

    l0 = frac(0.606433 + 1336.855225 * t)  # mean longitude of moon
    l = 2 * pi * frac(0.374897 + 1325.552410 * t)  # mean anomaly of Moon
    ls = 2 * pi * frac(0.993133 + 99.997361 * t)  # mean anomaly of Sun
    d = 2 * pi * frac(0.827361 + 1236.853086 * t)  # difference in longitude of moon and sun
    f = 2 * pi * frac(0.259086 + 1342.227825 * t)  # mean argument of latitude

    # corrections to mean longitude in arcsec
    dl = 22640 * sin(l)
    dl += -4586 * sin(l - 2 * d)
    dl += +2370 * sin(2 * d)
    dl += +769 * sin(2 * l)
    dl += -668 * sin(ls)
    dl += -412 * sin(2 * f)
    dl += -212 * sin(2 * l - 2 * d)
    dl += -206 * sin(l + ls - 2 * d)
    dl += +192 * sin(l + 2 * d)
    dl += -165 * sin(ls - 2 * d)
    dl += -125 * sin(d)
    dl += -110 * sin(l + ls)
    dl += +148 * sin(l - ls)
    dl += -55 * sin(2 * f - 2 * d)

    # simplified form of the latitude terms
    s = f + (dl + 412 * sin(2 * f) + 541 * sin(ls)) / arc
    h = f - 2 * d
    n = -526 * sin(h)
    n += +44 * sin(l + h)
    n += -31 * sin(-l + h)
    n += -23 * sin(ls + h)
    n += +11 * sin(-ls + h)
    n += -25 * sin(-2 * l + f)
    n += +21 * sin(-l + f)

    # ecliptic long and lat of Moon in rads
    l_moon = 2 * pi * frac(l0 + dl / 1296000)
    b_moon = (18520.0 * sin(s) + n) / arc

    # equatorial coord conversion - note fixed obliquity
    cb = cos(b_moon)
    x = cb * cos(l_moon)
    v = cb * sin(l_moon)
    w = sin(b_moon)
    y = coseps * v - sineps * w
    z = sineps * v + coseps * w
    # rho = sqrt(1.0 - z * z)
    # dec = (360.0 / 2 * pi) * atan(z / rho)
    # ra = ((48.0 / (2 * pi)) * atan(y / (x + rho))) % 24
    
    
    #########################
    # מרחק גאוצנטרי של הירח (בערך בקילומטרים)
    # הערכה פשוטה לפי הפרש באורך בין ירח ופריהליון
    # אפשרות מדויקת יותר: הוקטור עצמו *יחסי*, נחשב Δ לפי נורמה
    #distance_km = 385000.56 + dl  # קירוב, אפשר לשפר
    distance_km = sqrt(x**2 + y**2 + z**2) * 384400  # כפל ביחידה ממוצעת בק"מ
    #######################
    
    return x, y, z, distance_km  ################## המרחק זו תוספת שלי


###########################################################3

# פונקצייה מאוד חשובה מאת צ'אט גיפיטי להמרת קואורטינדות של גרם שמיימי מגאוצנטרי לטופוצנטרי כלומר ממיקום הצופה
def topocentric(x_geocentric, y_geocentric, z_geocentric, distance_km, lat_deg, lon_deg, lst_deg):
    
    # קבועים
    Re_km = 6378.137  # רדיוס כדור הארץ באקווטור בק"מ

    # המר לחישובים
    φ = radians(lat_deg)
    H = radians(lst_deg)  # זמן סידריאלי מקומי – זווית השעה

    # מיקום הצופה (יחידות רדיוס ארצי)
    x_obs = cos(φ) * cos(H)
    y_obs = cos(φ) * sin(H)
    z_obs = sin(φ)

    # המרחק לגוף השמימי ביחידות רדיוס ארצי
    rho = distance_km / Re_km

    # הפרש וקטור בין צופה לגוף השמימי
    xt = x_geocentric * rho - x_obs
    yt = y_geocentric * rho - y_obs
    zt = z_geocentric * rho - z_obs

    # נורמליזציה לווקטור יחידה
    r = sqrt(xt**2 + yt**2 + zt**2)
    
    x_topocentric = xt / r
    y_topocentric = yt / r
    z_topocentric = zt / r
    
    return x_topocentric, y_topocentric, z_topocentric

#############################################################


class RiSet:
    verbose = True
    # Riset.mtime() returns machine time as an int. The class variable tim is for
    # test purposes only and allows the hardware clock to be overridden
    tim = None
    ## For sunrise and sunset time search: What is the height of the sun above the horizon at sunrise or sunset.
    ##The normal default is -0.833 for normal sunrise and sunset
    sinho_sun_riset = -0.833 ## 

    @classmethod
    def mtime(cls):
        return round(time.time()) if cls.tim is None else cls.tim

    @classmethod
    def set_time(cls, t):  # Given time from Unix epoch set time
        if time.gmtime(0)[0] == 2000:  # Machine epoch
            t -= 10957 * 86400
        cls.tim = t

    def __init__(self, lat=LAT, long=LONG, lto=0, tl=None, dst=lambda x: x):  # Local defaults
        self.sglat = sin(radians(lat))
        self.cglat = cos(radians(lat))
        self.lat = lat ###########
        self.long = long
        self.check_lto(lto)  # -15 < lto < 15
        self.lto = round(lto * 3600)  # Localtime offset in secs
        self.tlight = sin(radians(tl)) if tl is not None else tl
        self.dst = dst
        self.mjd = None  # Current integer MJD
        # Times in integer secs from midnight on current day (in machine time adjusted for DST)
        # [sunrise, sunset, moonrise, moonset, cvend, cvstart]
        self._times = [None] * 6
        self.set_day()  # Initialise to today's date
        if RiSet.verbose:
            t = time.localtime()
            print(f"Machine time: {t[2]:02}/{t[1]:02}/{t[0]:4} {t[3]:02}:{t[4]:02}:{t[5]:02}")
            RiSet.verbose = False

    # ***** API start *****
    # Examine Julian dates either side of current one to cope with localtime.
    # 707μs on RP2040 at standard clock and with local time == UTC
    def set_day(self, day: int = 0):
        mjd = get_mjd(day)
        if self.mjd is None or self.mjd != mjd:
            spd = 86400  # Secs per day
            # ._t0 is time of midnight (local time) in secs since MicroPython epoch
            # time.time() assumes MicroPython clock is set to geographic local time
            self._t0 = ((self.mtime() + day * spd) // spd) * spd
            self.update(mjd)  # Recalculate rise and set times
        return self  # Allow r.set_day().sunrise()

    # variants: 0 secs since 00:00:00 localtime. 1 secs since MicroPython epoch
    # (relies on system being set to localtime). 2 human-readable text.
    def sunrise(self, variant: int = 0):
        return self._format(self._times[0], variant)

    def sunset(self, variant: int = 0):
        return self._format(self._times[1], variant)

    def moonrise(self, variant: int = 0):
        return self._format(self._times[2], variant)

    def moonset(self, variant: int = 0):
        return self._format(self._times[3], variant)

    def tstart(self, variant: int = 0):
        return self._format(self._times[4], variant)

    def tend(self, variant: int = 0):
        return self._format(self._times[5], variant)

    def set_lto(self, t):  # Update the offset from UTC
        self.check_lto(t)  # No need to recalc beause date is unchanged
        self.lto = round(t * 3600)  # Localtime offset in secs

    def has_risen(self, sun: bool):
        return self.has_x(True, sun)

    def has_set(self, sun: bool):
        return self.has_x(False, sun)

    # Return current state of sun or moon. The moon has a special case where it
    # rises and sets in a 24 hour period. If its state is queried after both these
    # events or before either has occurred, the current state depends on the order
    # in which they occurred (the sun always sets afer it rises).
    # The case is (.has_risen(False) and .has_set(False)) and if it occurs then
    # .moonrise() and .moonset() must return valid times (not None).
    def is_up(self, sun: bool):
        hr = self.has_risen(sun)
        hs = self.has_set(sun)
        rt = self.sunrise() if sun else self.moonrise()
        st = self.sunset() if sun else self.moonset()
        if rt is None and st is None:  # No event in 24hr period.
            return self.above_horizon(sun)
        # Handle special case: moon has already risen and set or moon has neither
        # risen nor set, yet there is a rise and set event in the day
        if not (hr ^ hs):
            if not ((rt is None) or (st is None)):
                return rt > st
        if not (hr or hs):  # No event has yet occurred
            return rt is None

        return hr and not hs  # Default case: up if it's risen but not set

    # ***** API end *****

    # Generic has_risen/has_set function
    def has_x(self, risen: bool, sun: bool):
        if risen:
            st = self.sunrise(1) if sun else self.moonrise(1)  # DST- adjusted machine time
        else:
            st = self.sunset(1) if sun else self.moonset(1)
        if st is not None:
            return st < self.dst(self.mtime())  # Machine time
        return False

    def above_horizon(self, sun: bool):
        now = self.mtime() + self.lto  # UTC
        tutc = (now % 86400) / 3600  # Time as UTC hour of day (float)
        return self.sin_alt(tutc, sun) > 0  # Object is above horizon

    # Re-calculate rise and set times
    def update(self, mjd):
        for x in range(len(self._times)):
            self._times[x] = None  # Assume failure
        days = (1, 2) if self.lto < 0 else (1,) if self.lto == 0 else (0, 1)
        tr = None  # Assume no twilight calculations
        ts = None
        for day in days:
            self.mjd = mjd + day - 1
            sr, ss = self.rise_set(True, False)  # Sun
            # Twilight: only calculate if required
            if self.tlight is not None:
                tr, ts = self.rise_set(True, True)
            mr, ms = self.rise_set(False, False)  # Moon
            # Adjust for local time and DST. Store in ._times if value is in
            # 24-hour local time window
            self.adjust((sr, ss, mr, ms, tr, ts), day)
        self.mjd = mjd

    def adjust(self, times, day):
        for idx, n in enumerate(times):
            if n is not None:
                n += self.lto + (day - 1) * 86400
                n = self.dst(n)  # Adjust for DST on day of n
                h = n // 3600
                if 0 <= h < 24:
                    self._times[idx] = n

    def _format(self, n, variant):
        if (n is not None) and (variant & 4):  # Machine clock set to UTC
            variant &= 0x03
            n = self.dst(n + self._t0) - self._t0
        if variant == 0:  # Default: secs since Midnight (local time)
            return n
        elif variant == 1:  # Secs since epoch of MicroPython platform
            return None if n is None else n + self._t0
        # variant == 2
        if n is None:
            return "--:--:--"
        else:
            hr, tmp = divmod(n, 3600)
            mi, sec = divmod(tmp, 60)
            return f"{hr:02d}:{mi:02d}:{sec:02d}"

    def check_lto(self, t):
        if not -15 < t < 15:
            raise ValueError("Invalid local time offset.")

    # See https://github.com/orgs/micropython/discussions/13075
    def lstt(self, t, h):
        # Takes the mjd and the longitude (west negative) and then returns
        # the local sidereal time in degrees. Im using Meeus formula 11.4
        # instead of messing about with UTo and so on
        # modified to use the pre-computed 't' value from sin_alt
        d = t * 36525
        df = frac(0.5 + h / 24)
        c1 = 360
        c2 = 0.98564736629
        dsum = c1 * df + c2 * d  # dsum is still ~ 9000 on average, losing precision
        lst = 280.46061837 + dsum + t * t * (0.000387933 - t / 38710000)
        return lst

    def sin_alt(self, hour, sun):
        # Returns the sine of the altitude of the object (moon or sun)
        # at an hour relative to the current date (mjd)
        func = minisun if sun else minimoon
        mjd = (self.mjd - 51544.5) + hour / 24.0
        # mjd = self.mjd + hour / 24.0
        t = mjd / 36525.0
        x, y, z, distance_km = func(t) ################# הוספתי את distance_km
        tl = self.lstt(t, hour) + self.long  # Local mean sidereal time adjusted for logitude
        return self.sglat * z + self.cglat * (x * cos(radians(tl)) + y * sin(radians(tl)))
    


    ######################################################################################################3
    #######################################################################################################


    # פונקצייה שבניתי על בסיס הפונקצייה הקודמת ביחד עם בינה מלאכותית ודוגמאות מספריית אסטרל
    # הפונקצייה מחזירה גובה השמש או הירח במעלות, אזימוט מהצפון במעלות, עלייה ישרה בשעות (עשרוני-ללא דקות ושניות), דקלינציה במעלות
    # בתוך הפונקצייה מחושב גם זמן כוכבים מקומי וזווית השעה של השמש או הירח אך הם לא מוחזרים
    def alt_az_ra_dec(self, hour, sun=True):
        """
        מחזירה את גובה השמש (Alt) ואת האזימוט שלה (Az) במעלות.
        """
        func = minisun if sun else minimoon
        mjd = (self.mjd - 51544.5) + hour / 24.0
        t = mjd / 36525.0
        tl = self.lstt(t, hour) + self.long  # זמן כוכבי מקומי במעלות
        
        # זה מחזיר את מיקום הגוף השמיימי מנקודת מבט גיאוצנטרית כלומר ממרכז כדור הארץ
        xg, yg, zg, distance_km = func(t)
        
        # זה ממיר לטופוצנטרי לצורך חישוב עתידי מדוייק יותר של גובה
        xt, yt, zt = topocentric(xg, yg, zg, distance_km, self.lat, self.long, tl)

        # חישוב גובה מהאופק על מיקום טופוצנטרי של הצופה
        sin_alt = self.sglat * zt + self.cglat * (xt * cos(radians(tl)) + yt * sin(radians(tl)))
        alt = degrees(asin(sin_alt))  # גובה השמש במעלות
        
        # חישוב נטייה כלומר דקלינציה וזה דווקא על מיקום גיאוצנטרי ממרכז כדור הארץ
        rho = sqrt(xg * xg + yg * yg)  # היטל של הווקטור על מישור XY
        dec = degrees(atan2(zg, rho))  # חישוב נטייה במעלות
        
        # חישוב עלייה ישרה (RA) וזה דווקא גיאוצנטרי ממרכז כדור הארץ
        # בתאריך 14.3.25 בשעה 21:11:00 utc+2 במודיעין עילית קיבלתי שגיאת חלוקה באפס
        # השגיאה מגיעה מ x + rho ולהלן הפתרון שצ'אט גיפיטי הציע
        ra_base = xg + rho
        epsilon = 1e-9  # ערך קטן מאוד שקרוב מאוד ל 0.00000 אבל הוא לא אפס מוחלט וזאת כדי למנוע שגיאת חלוקה באפס 
        if abs(ra_base) < epsilon:  
            ra_base = epsilon  # לשמור על דיוק גבוה אבל למנוע בעיה
        ra = ((48.0 / (2 * pi)) * atan(yg / ra_base)) % 24 # עלייה ישרה בשעות כשבר עשרוני

        # חישוב האזימוט (Az) מתוך העלייה הישרה
        hourangle = radians(tl) - radians(ra * 15)  # זמן הכוכבים המקומי פחות העלייה הישרה של הכוכב זה זוית השעה שלו (ra * 15 מחזיר למעלות)
        hourangle_hours = (degrees(hourangle) % 360)  / 15.0 # זווית השעה בשעות כשבר עשרוני
        
        sh = sin(hourangle)
        ch = cos(hourangle)
        sd = sin(radians(dec))
        cd = cos(radians(dec))
        sl = self.sglat # ==sin(radians(lat))
        cl = self.cglat # ==cos(radians(lat))

        x_az = -ch * cd * sl + sd * cl
        y_az = -sh * cd
        az = degrees(atan2(y_az, x_az)) % 360  # אזימוט במעלות
        
        return alt, az, ra, dec

    

######################################################################################################3
#######################################################################################################

    # Calculate rise and set times of sun or moon for the current MJD. Times are
    # relative to that 24 hour period.
    def rise_set(self, sun, tl):
        t_rise = None  # Rise and set times in secs from midnight
        t_set = None
        if tl:
            sinho = -self.tlight
        else:
            sinho = sin(radians(RiSet.sinho_sun_riset)) if sun else sin(radians(8 / 60)) ## sinho_sun_riset
        # moonrise taken as centre of moon at +8 arcmin
        # sunset upper limb simple refraction
        # The loop finds the sin(alt) for sets of three consecutive
        # hours, and then tests for a single zero crossing in the interval
        # or for two zero crossings in an interval for for a grazing event
        yp = self.sin_alt(0, sun) - sinho
        for hour in range(1, 24, 2):
            ym = yp
            yz = self.sin_alt(hour, sun) - sinho
            yp = self.sin_alt(hour + 1, sun) - sinho
            nz, z1, z2, ye = quad(ym, yz, yp)  # Find horizon crossings
            if nz == 1:  # One crossing found
                if ym < 0.0:
                    t_rise = 3600 * (hour + z1)
                else:
                    t_set = 3600 * (hour + z1)
            # case where two events are found in this interval
            # (rare but whole reason we are not using simple iteration)
            elif nz == 2:
                if ye < 0.0:
                    t_rise = 3600 * (hour + z2)
                    t_set = 3600 * (hour + z1)
                else:
                    t_rise = 3600 * (hour + z1)
                    t_set = 3600 * (hour + z2)

            if t_rise is not None and t_set is not None:
                break  # All done
        return to_int(t_rise), to_int(t_set)  # Convert to int preserving None values


############################################################################################################################
                                # 2. moonphase_sgb
############################################################################################################################
                                

# moonphase.py Calculate lunar phases

# Source Tech\ Notes/Astronomy/astro_references/moontool.c
# The information for this was drawn from public domain sources including C code
# written by John Walker and Ron Hitchens in 1987-88 and released with the "licence"
# Do what thou wilt shall be the whole of the law".

# Uses Python arbitrary length integers to maintain accuracy on platforms with
# 32-bit floating point.
# Copyright (c) Peter Hinch 2023  Released under the MIT license.

## Changes made by Simcha Gershon Bohrer marked with ##


# Exports calc_phases()

from math import radians, sin, cos, floor
import time
import array

SYNMONTH = 29.53058868  # Synodic month (new Moon to new Moon)


# MEANPHASE  --  Calculates time of the mean new Moon for a given base date.
# This argument K to this function is the precomputed synodic month index, given by:
# K = (year - 1900) * 12.3685
# where year is expressed as a year and fractional year.
# sdate is days from 1900 January 0.5. Returns days from 1900 January 0.5


def meanphase(sdate: float, k: int) -> float:
    # Time in Julian centuries from 1900 January 0.5
    t = sdate / 36525
    t2 = t * t  # Square for frequent use
    t3 = t2 * t  # Cube for frequent use
    nt1 = 0.75933 + SYNMONTH * k + 0.0001178 * t2 - 0.000000155 * t3
    return nt1 + 0.00033 * sin(radians(166.56 + 132.87 * t - 0.009173 * t2))


# TRUEPHASE  --  Given a K value used to determine the mean phase of the new moon,
# and a phase no.  (0..3), return the true, corrected phase time
# as integer Julian seconds.


def truephase(k: int, phi: int) -> int:
    k += (0, 0.25, 0.5, 0.75)[phi]  # Add phase to new moon time
    t = k / 1236.85  # Time in Julian centuries from 1900 January 0.5
    t2 = t * t  # Square for frequent use
    t3 = t2 * t  # Cube for frequent use
    # Sun's mean anomaly
    m = 359.2242 + 29.10535608 * k - 0.0000333 * t2 - 0.00000347 * t3
    # Moon's mean anomaly
    mprime = 306.0253 + 385.81691806 * k + 0.0107306 * t2 + 0.00001236 * t3
    # Moon's argument of latitude
    f = 21.2964 + 390.67050646 * k - 0.0016528 * t2 - 0.00000239 * t3
    if phi in (0, 2):  # Corrections for New and Full Moon
        pt = (0.1734 - 0.000393 * t) * sin(radians(m))
        pt += 0.0021 * sin(radians(2 * m))
        pt -= 0.4068 * sin(radians(mprime))
        pt += 0.0161 * sin(radians(2 * mprime))
        pt -= 0.0004 * sin(radians(3 * mprime))
        pt += 0.0104 * sin(radians(2 * f))
        pt -= 0.0051 * sin(radians(m + mprime))
        pt -= 0.0074 * sin(radians(m - mprime))
        pt += 0.0004 * sin(radians(2 * f + m))
        pt -= 0.0004 * sin(radians(2 * f - m))
        pt -= 0.0006 * sin(radians(2 * f + mprime))
        pt += 0.0010 * sin(radians(2 * f - mprime))
        pt += 0.0005 * sin(radians(m + 2 * mprime))
    else:  # First or last quarter
        pt = (0.1721 - 0.0004 * t) * sin(radians(m))
        pt += 0.0021 * sin(radians(2 * m))
        pt -= 0.6280 * sin(radians(mprime))
        pt += 0.0089 * sin(radians(2 * mprime))
        pt -= 0.0004 * sin(radians(3 * mprime))
        pt += 0.0079 * sin(radians(2 * f))
        pt -= 0.0119 * sin(radians(m + mprime))
        pt -= 0.0047 * sin(radians(m - mprime))
        pt += 0.0003 * sin(radians(2 * f + m))
        pt -= 0.0004 * sin(radians(2 * f - m))
        pt -= 0.0006 * sin(radians(2 * f + mprime))
        pt += 0.0021 * sin(radians(2 * f - mprime))
        pt += 0.0003 * sin(radians(m + 2 * mprime))
        pt += 0.0004 * sin(radians(m - 2 * mprime))
        pt -= 0.0003 * sin(radians(2 * m + mprime))
        if phi < 2:  # First quarter correction
            pt += 0.0028 - 0.0004 * cos(radians(m)) + 0.0003 * cos(radians(mprime))
        else:  # Last quarter correction
            pt += -0.0028 + 0.0004 * cos(radians(m)) - 0.0003 * cos(radians(mprime))
    pt = round(pt * 86400)  # Integer seconds from here
    pt += round(2_953_058_868 * 864 * k) // 1000_000  # round(SYNMONTH * k * 86400)
    qq = 0.0001178 * t2 - 0.000000155 * t3
    qq += 0.00033 * sin(radians(166.56 + 132.87 * t - 0.009173 * t2))
    pt += round(qq * 86400)  # qq amounts to 2s
    return pt + 208_657_793_606


def dt_to_text(tim):  # Convert a time to text
    t = time.localtime(tim)
    return f"{t[2]:02}/{t[1]:02}/{t[0]:4} {t[3]:02}:{t[4]:02}:{t[5]:02}"


class MoonPhase:
    verbose = True

    #######################################################################################
    # כל זה הוספתי כמו ב sun_moon כדי שיהיה אפשר להגדיר זמן אחר במקום זמן המכונה הפנימי
    # MoonPhase.mtime() returns machine time as an int. The class variable tim is for
    # test purposes only and allows the hardware clock to be overridden
    tim = None 

    @classmethod
    def mtime(cls):
        return round(time.time()) if cls.tim is None else cls.tim

    
    @classmethod
    def set_time(cls, t):  # Given time from Unix epoch set time
        if time.gmtime(0)[0] == 2000:  # Machine epoch
            t -= 10957 * 86400
        cls.tim = t
    
    ##################################################################################
    
    def __init__(self, lto: float = 0, dst=lambda x: x):
        self.lto_s = self._check_lto(lto)  # -15 < lto < 15
        # local time = UTC + lto .lto_s = offset in secs
        self.dst = dst
        # Datetimes in secs since hardware epoch based on UTC
        # With epoch 1970 this could need long ints.
        self.phases = array.array("q", (0,) * 5)
        # Calculate Julian date of machine epoch
        # Multiply by 100 to avoid fraction
        jepoch = 244058750  # Julian date of Unix epoch (1st Jan 1970) * 100
        if time.gmtime(0)[0] == 2000:  # Machine epoch
            jepoch += 1095700
        jepoch *= 864  # Seconds from epoch
        self.jepoch = jepoch
        self.secs = 0  # Time of calling .set_day in secs UTC
        self.set_day()  # Populate array and .secs
        if MoonPhase.verbose:
            print(f"Machine time: {dt_to_text(time.time())}")
            MoonPhase.verbose = False

    # Take offset in days from today, return time of last midnight in secs from machine epoch
    # Take time of last midnight machine time in secs since machine epoch. Add a
    # passed offset in days. Convert to UTC using LTO. The returned value is as
    # if the hardware clock were running UTC.
    def _midnight(self, doff: float = 0):  # Midnight last night + days offset (UTC)
        tl = round((self.mtime() // 86400 + doff) * 86400)  # Target in local time     ## self.mtime() במקום time.time()
        return tl - self.lto_s

    def set_lto(self, t: float):  # Update the offset from UTC
        self.lto_s = self._check_lto(t)  # Localtime offset in secs

    def set_day(self, doff: float = 0):
        self.secs = round(self.mtime() + doff * 86400 - self.lto_s)   ## self.mtime() במקום time.time()
        start = self._midnight(doff)  # Phases are calculated around this time (UTC)
        self._populate(start)  # Immediate return if .phases already OK

    def datum(self, text: bool = True):
        t = self.secs + self.lto_s
        return dt_to_text(t) if text else t

    def quarter(self, q: int, text: bool = True):
        if not 0 <= q <= 4:
            raise ValueError("Quarter nos must be from 0 to 4.")
        tutc = self.phases[q]  # Time of phase in secs UTC
        # Adjust time: t is machine time in secs since machine epoch
        t = self.dst(tutc + self.lto_s)  # UTC secs from hardware epoch -> local time
        return dt_to_text(t) if text else t  # Secs since machine epoch

    # Return moon  phase as 0.0 <= n < 1.0 by defaut for current datetime.
    def phase(self) -> float:  # doff: days offset with optional fraction
        t = self.secs  # As set by .set_day()
        if not (self.phases[0] <= t <= self.phases[4]):  # set_day was not called
            self.set_day()  # Assume today
        prev = self.phases[0]
        for n, phi in enumerate(self.phases):
            if phi > t:
                break  # phi is upcoming phase time
            prev = phi  # Last phase before now
        if prev == phi:  # Day is day of new moon: use synodic month/4
            r = (t - prev) * 0.25 / 637860.715488
            if r < 0:
                r = 1 - r
        else:
            r = (n - 1) * 0.25 + (t - prev) * 0.25 / (phi - prev)
        return min(r, 0.999999)  # Rare pathological results where r slightly > 1.0

    def _next_lunation(self):  # Use approx time of next full moon to advance
        self._populate(round(self.phases[2] + SYNMONTH * 86400))

    # toff: days offset with optional fraction
    def nextphase(self, text: bool = True):
        n = 0
        lun = 0  # Skip historic quarters
        while True:
            yield n, lun, self.quarter(n, text)
            n += 1
            n %= 4
            if n == 0:
                self._next_lunation()
                lun += 1

    def _check_lto(self, lto: float) -> int:
        if not -15 < lto < 15:
            raise ValueError("Invalid local time offset.")
        return round(lto * 3600)

    # Populate the phase array. Fast return if phases are alrady correct.
    # Find time of phases of the moon which surround the passed datetime.
    # Five phases are found, starting and ending with the new moons which bound
    # the specified lunation.
    # Passed time, and the result in .phases, are seconds since hardware epoch
    # adjusted for UTC: i.e. as if the RTC were running UTC rather than local time.
    def _populate(self, t: int):
        if self.phases[0] < t < self.phases[4]:
            return  # Nothing to do
        # Return days since Jan 0.5 1900 as a float. Returns same value on 32 and 64 bits
        def jd1900(t: int) -> float:
            y, m, mday = time.localtime(t)[:3]
            if m <= 2:
                m += 12
                y -= 1
            b = round(y / 400 - y / 100 + y / 4)
            mjm = 365 * y - 679004 + b + int(30.6001 * (m + 1)) + mday
            return mjm - 15019.5

        sdate: float = jd1900(t)  # Days since 1900 January 0.5
        adate: float = sdate - 45
        yy, mm, dd = time.localtime(t)[:3]
        k1: int = floor((yy + ((mm - 1) * (1.0 / 12.0)) - 1900) * 12.3685)  # 365.25/SYNMONTH
        adate = meanphase(adate, k1)  # Find new moon well before current date
        nt1: float = adate
        while True:
            adate += SYNMONTH  # For each lunar month
            k2: int = k1 + 1
            nt2: float = meanphase(adate, k2)
            if nt1 <= sdate and nt2 > sdate:
                break
            nt1 = nt2
            k1 = k2
        # k is integer days since start of 1900, being the lunation number
        # 1533, 1534 on both platforms.
        for n, k in enumerate((k1, k1, k1, k1, k2)):
            phi: int = truephase(k, n % 4)  # Args lunation no., phase no. 0..3
            self.phases[n] = phi - self.jepoch  # Julian datetime to secs since hardware epoch
            # Datetimes in secs since hardware epoch based on UTC

###################################################################################################################################
                            # 3. mpy_heb_date
###################################################################################################################################

import time


###################################################################
       # גמטריה מספריית פיילוח לצורך שנה עברית באותיות
###################################################################

# https://github.com/simlist/pyluach/blob/master/src/pyluach/gematria.py

_GEMATRIOS = {
    1: 'א',
    2: 'ב',
    3: 'ג',
    4: 'ד',
    5: 'ה',
    6: 'ו',
    7: 'ז',
    8: 'ח',
    9: 'ט',
    10: 'י',
    20: 'כ',
    30: 'ל',
    40: 'מ',
    50: 'נ',
    60: 'ס',
    70: 'ע',
    80: 'פ',
    90: 'צ',
    100: 'ק',
    200: 'ר',
    300: 'ש',
    400: 'ת'
}


def _stringify_gematria(letters):
    """Insert geresh or gershayim symbols into gematria."""
    length = len(letters)
    if length > 1:
        return f'{letters[:-1]}״{letters[-1]}'
    if length == 1:
        return f'{letters}׳'
    return ''


def _get_letters(num):
    """Convert numbers under 1,000 into raw letters."""
    ones = num % 10
    tens = num % 100 - ones
    hundreds = num % 1000 - tens - ones
    four_hundreds = ''.join(['ת' for i in range(hundreds // 400)])
    ones = _GEMATRIOS.get(ones, '')
    tens = _GEMATRIOS.get(tens, '')
    hundreds = _GEMATRIOS.get(hundreds % 400, '')
    letters = f'{four_hundreds}{hundreds}{tens}{ones}'
    return letters.replace('יה', 'טו').replace('יו', 'טז')


def _num_to_str(num, thousands=False, withgershayim=True):
    """Return gematria string for number.

    Parameters
    ----------
    num : int
        The number to get the Hebrew letter representation
    thousands : bool, optional
        True if the hebrew returned should include a letter for the
        thousands place ie. 'ה׳' for five thousand.

    Returns
    -------
    str
        The Hebrew representation of the number.
    """
    letters = _get_letters(num)
    if withgershayim:
        letters = _stringify_gematria(letters)
    if thousands:
        thousand = _get_letters(num // 1000)
        if withgershayim:
            thousand = ''.join([thousand, '׳'])
        letters = ''.join([thousand, letters])
    return letters


###################################################################
                    # עד כאן גמטריה מספריית פיילוח
###################################################################

# פונקצייה עבור mktime שמחזירה טאפל באורך 9 מקומות שזה מתאים לפייתון רגיל או מיקרופייתון
# המיקום התשיעי אומר אם זה שעון קיץ או חורף לפי אפס או אחד. ומינוס אחד אומר שהמחשב יחליט לבד אם זה שעון קיץ או חורף
def get_mktime_date_only(year, month, day):
    return time.mktime((year, month, day, 0, 0, 0, 0, 0, -1))
    

# מוציא את מספר היום בשבוע הנורמלי לפי סדר מתוך שעון המכשיר שמוגדר RTC
def get_normal_weekday(rtc_weekday):
    weekday_dict = {6:1,0:2,1:3,2:4,3:5,4:6,5:7}
    return weekday_dict.get(rtc_weekday)


def get_holiday_name(heb_day_int, heb_month_int, is_leap_year):
    """ מקבלת יום, חודש והאם השנה מעוברת, ומחזירה את שם החג אם מדובר בחג, אחרת מחזירה False """
    HOLIDAYS = {
        (1, 1): "ראש השנה",
        (10, 1): "יום כיפור",
        (15, 1): "ראשון של סוכות",
        (22, 1): "שמיני עצרת",
        (15, 8 if is_leap_year else 7): "ראשון של פסח",
        (21, 8 if is_leap_year else 7): "שביעי של פסח",
        (6, 10 if is_leap_year else 9): "שבועות"  
    }
    return HOLIDAYS.get((heb_day_int, heb_month_int), False)

def get_lite_holiday_name(heb_day_int, heb_month_int, is_leap_year, is_kislev_29):
    """ מקבלת יום, חודש והאם השנה מעוברת, ומחזירה את שם החג הקל כלומר מדרבנן אם מדובר בחג קל, אחרת מחזירה False """
    
    #  שימו לב שהצומות יכולים להתאחר ליום ראשון אם הם חלים בשבת
    # הם גם יכולים להקדים ליום חמישי אם הם חלים ביום שישי למעט עשרה בטבת שלא מקדימים אותו 
    # אבל הקדמת ואיחור הצומות לא מטופל כאן כרגע
    
    LITE_HOLIDAYS = {

        (16, 1): "א דחול המועד סוכות",
        (17, 1): "ב דחול המועד סוכות",
        (18, 1): "ג דחול המועד סוכות",
        (19, 1): "ד דחול המועד סוכות",
        (20, 1): "ה דחול המועד סוכות",
        (21, 1): "הושענה רבא",
        (25, 3): "נר ראשון של חנוכה",
        (26, 3): "נר שני של חנוכה",
        (27, 3): "נר שלישי של חנוכה",
        (28, 3): "נר רביעי של חנוכה",
        (29, 3): "נר חמישי של חנוכה",
        (1, 4) if is_kislev_29 else (30, 3): "נר שישי של חנוכה",
        (2, 4) if is_kislev_29 else (1, 4): "נר שביעי של חנוכה",
        (3, 4) if is_kislev_29 else (2, 4): "נר שמיני של חנוכה",
        (14, 7) if is_leap_year else (14, 6): "פורים דפרזים",
        (15, 7) if is_leap_year else (15, 6): "פורים דמוקפין",
        (16, 8) if is_leap_year else (16, 7): "א דחול המועד פסח",
        (17, 8) if is_leap_year else (17, 7): "ב דחול המועד פסח",
        (18, 8) if is_leap_year else (18, 7): "ג דחול המועד פסח",
        (19, 8) if is_leap_year else (19, 7): "ד דחול המועד פסח",
        (20, 8) if is_leap_year else (20, 7): "ה דחול המועד פסח",
        (3, 1): "צום גדליה",
        (10, 4): "צום עשרה בטבת", 
        (17, 11) if is_leap_year else (17, 10): "צום שבעה עשר בתמוז", 
        (9, 12) if is_leap_year else (9, 11): "צום תשעה באב",    
    }
    return LITE_HOLIDAYS.get((heb_day_int, heb_month_int), False)


# מילון לשמות החודשים בעברית
def heb_month_names(number, is_leap=False):
    d={
        1:"תשרי",
        2:"מרחשוון",
        3:"כסלו",
        4:"טבת",
        5:"שבט",
        6:"אדר" if not is_leap else "אדר-א",
        7:"ניסן" if not is_leap else "אדר-ב",
        8:"אייר" if not is_leap else "ניסן",
        9:"סיוון" if not is_leap else "אייר",
        10:"תמוז" if not is_leap else "סיוון",
        11:"אב" if not is_leap else "תמוז",
        12:"אלול" if not is_leap else "אב",
        13:"" if not is_leap else "אלול",}
    return d.get(number)

# מילון לשמות הימים בחודש בעברית
def heb_month_day_names(number):
    d={
        1:"א",
        2:"ב",
        3:"ג",
        4:"ד",
        5:"ה",
        6:"ו",
        7:"ז",
        8:"ח",
        9:"ט",
        10:"י",
        11:"יא",
        12:"יב",
        13:"יג",
        14:"יד",
        15:"טו",
        16:"טז",
        17:"יז",
        18:"יח",
        19:"יט",
        20:"כ",
        21:"כא",
        22:"כב",
        23:"כג",
        24:"כד",
        25:"כה",
        26:"כו",
        27:"כז",
        28:"כח",
        29:"כט",
        30:"ל",}
    return d.get(number)

# מילון לשמות הימים בשבוע בעברית
def heb_weekday_names(number):
    d={
        1:"ראשון",
        2:"שני",
        3:"שלישי",
        4:"רביעי",
        5:"חמישי",
        6:"שישי",
        7:"שבת",}
    return d.get(number)


# מילון למבני השנים האפשריים בלוח העברי לפי מספר ימי השנה נותן את מספר הימים שיש בכל חודש
def get_year_structure(year_length):
    
    # מבני השנים האפשריים
    structures = {
        353: [30, 29, 29, 29, 30, 29, 30, 29, 30, 29, 30, 29],
        354: [30, 29, 30, 29, 30, 29, 30, 29, 30, 29, 30, 29],
        355: [30, 30, 30, 29, 30, 29, 30, 29, 30, 29, 30, 29],
        383: [30, 29, 29, 29, 30, 30, 29, 30, 29, 30, 29, 30, 29],
        384: [30, 29, 30, 29, 30, 30, 29, 30, 29, 30, 29, 30, 29],
        385: [30, 30, 30, 29, 30, 30, 29, 30, 29, 30, 29, 30, 29]
    }
    return structures.get(year_length)

# פונקצייה נורא חשובה שמקבלת קלט של תאריך עברי שממנו רוצים להזיז ימים וקלט של כמה ימים רוצים להזיז וקלט מהו אורך השנה העברית
# ואז היא אומרת לאיזה תאריך הגענו. היא נבנתה רק על ידי צאט גיפיטי על בסיס נתונים שנתתי לו
def move_heb_date(start_day, start_month, year_length, days_to_move):
    # קבלת מבנה השנה
    year_structure = get_year_structure(year_length)
    if not year_structure:
        raise ValueError("אורך השנה לא תקין")

    # האם השנה מעוברת
    is_leap = year_length in [383, 384, 385]

    # חישוב היום החדש
    current_day = start_day
    current_month = start_month

    # הזזה קדימה או אחורה
    while days_to_move != 0:
        days_in_month = year_structure[current_month - 1]
        if days_to_move > 0:  # הזזה קדימה
            remaining_days_in_month = days_in_month - current_day
            if days_to_move <= remaining_days_in_month:
                current_day += days_to_move
                days_to_move = 0
            else:
                days_to_move -= (remaining_days_in_month + 1)
                current_day = 1
                current_month += 1
                if current_month > len(year_structure):  # מעבר לשנה הבאה
                    if days_to_move == 0:  # בדיוק ביום האחרון
                        current_month -= 1
                        current_day = year_structure[current_month - 1]
                    else:
                        raise ValueError("החישוב חרג מגבולות השנה")
        else:  # הזזה אחורה
            if abs(days_to_move) < current_day:
                current_day += days_to_move
                days_to_move = 0
            else:
                days_to_move += current_day
                current_month -= 1
                if current_month < 1:  # מעבר לשנה קודמת
                    raise ValueError("החישוב חרג מגבולות השנה")
                current_day = year_structure[current_month - 1]
                
    return current_day, current_month




# פונקצייה שמחזירה את התאריך הגרגוריאני שבו יחול פסח בשנה נתונה או את התאריך הגרגוריאני שבו יחול ראש השנה שאחרי פסח של השנה הנתונה
# כברירת מחדל מקבל קלט של שנה לועזית אך יכול לקבל קלט של שנה עברית במספרים אם מגדירים זאת בקריאה לפונקצייה
def get_geus_rosh_hashana_greg(year, from_heb_year = False):

    if from_heb_year:
        A = year
        # הגדרת שנה לועזית המקבילה לשנה העברית שהוזנה
        B = A - 3760

    else:
        B = year
        A = B + 3760

    # אינני יודע מה מייצגות שתי ההגדרות הבאות 

    # איי קטנה נותן מספר בין 0 ל- 18 שממנו יודעים האם השנה העברית פשוטה או מעוברת. אם איי קטנה קטן מ-11 השנה היא פשוטה, ואם גדול מ-12 השנה היא מעוברת
    # בנוסף, ככל שאיי קטנה קרובה יותר למספר 18, זה אומר שפסח רחוק יותר מתקופת ניסן
    a = (12 * A + 17) % 19
    
    # נוסחה לקבל את מספר השנה במחזור השנים הפשוטות והמעוברות לפי איי קטנה
    # לדוגמא אם איי קטנה שווה 10 אז מספר השנה במחזור 19 השנים הוא 1
    shana_bemachzor19 = {10:1,3:2,15:3,8:4,1:5,13:6,6:7,18:8,11:9,4:10,16:11,9:12,2:13,14:14,7:15,0:16,12:17,5:18,17:19}.get(a)

    # בי קטנה מציינת האם השנה היוליאנית המקבילה היא פשוטה (365 יום) או כבושה (366 יום). אם אין שארית, השנה היא כבושה
    b = A % 4

    # נוסחת גאוס בשברים עשרוניים
    nuscha = 32.0440931611436 + (1.5542417966211826) * a + 0.25 * b - (0.0031777940220922675) * A 

    # נוסחת גאוס בשברים פשוטים
    #nuscha = 32 + 4343/98496 + (1 + 272953/492480) * a + 1/4 * b - (313/98496) * A

    # אם גדולה זה השלם של הנוסחה
    # ט"ו בניסן של השנה המבוקשת יחול ביום אם גדולה בחודש מרס
    M = int(nuscha)

    # אם קטנה היא השארית של הנוסחה, והיא חשובה לצורך הדחיות
    m = nuscha - int(nuscha)

    # סי הוא היום בשבוע שבו יחול פסח של השנה המבוקשת. אם סי שווה לאפס הכוונה ליום שבת 7
    c = (M + 3 * A + 5 * b + 5) % 7

    # מידע: דחיית מולד זקן מוכנסת כבר במספר 32 שבנוסחה הראשית

    # חישוב דחיית לא בד"ו פסח שהיא שיקוף של דחיית לא אד"ו ראש
    if c in (2,4,6):
        c = c + 1
        M = M + 1
    # חישוב השפעת דחיית גטר"ד בשנה פשוטה
    elif c == 1 and a > 6 and m >= 0.6329:
        c = c + 2
        M = M + 2
    # חישוב השפעת דחיית בטו תקפט בשנה פשוטה שהיא מוצאי מעוברת
    elif c == 0 and a > 11 and m >= 0.8977:
        c = c + 1
        M = M + 1
    else:
        c = c
        M = M

    # טיפול באם היום בשבוע של פסח יוצא אפס זה אומר יום 7 שזה שבת
    if c == 0:
        c = c + 7

    # אם אם גדולה קטן או שווה לשלושים ואחד פסח יהיה בחודש מרס
    if M <= 31:
        M = M
        chodesh_julyani_pesach = 3 
    # במצב הבא התאריך יהיה בחודש אפריל במקום בחודש מרס
    elif M > 31:
        M = M - 31
        chodesh_julyani_pesach = 4
        
        
    # מעבר ללוח הגרגוריאני
    # חודש מרס הוא תמיד 31 ימים

    if B >= 1582 and B < 1700:
        M = (M + 10) 
    elif B >= 1700 and B < 1800:
        M = (M + 11) 
    elif B >= 1800 and B < 1900:
        M = (M + 12) 
    elif B >= 1900 and B < 2100:
        M = (M + 13) 
    elif B >= 2100 and B < 2200:
        M = (M + 14) 
    elif B >= 2200 and B < 2300:
        M = (M + 15) 
    else:
        M = M

    # אם אם גדולה קטן או שווה לשלושים ואחד פסח יהיה בחודש מרס
    if M <= 31:
        M = M
        chodesh_gregoriani_pesach = chodesh_julyani_pesach

    # במצב הבא התאריך יהיה בחודש אפריל במקום בחודש מרס
    elif M > 31:
        M = M - 31
        chodesh_gregoriani_pesach = chodesh_julyani_pesach + 1

    pesach_greg_day = M
    pesach_greg_month = chodesh_gregoriani_pesach
    pesach_greg_year = B
    pesach_weekday = c
    
    # האם זו שנה עברית מעוברת
    heb_leap_year = shana_bemachzor19 in (3,6,8,11,14,17,19)
    
    #############################################################################################################
    # מציאת התאריך הלועזי של ראש השנה של השנה הבא לאחר הפסח ראו ספר שערים ללוח העברי עמוד 204
    next_rosh_hashana_greg_day = pesach_greg_day + 10
    if pesach_greg_month == 3:
        next_rosh_hashana_greg_month = 8
    elif pesach_greg_month == 4:
        next_rosh_hashana_greg_month = 9
        
    next_rosh_hashana_greg_year = pesach_greg_year
    
    if next_rosh_hashana_greg_day > 31 and pesach_greg_month == 3:
        next_rosh_hashana_greg_day = next_rosh_hashana_greg_day - 31
        next_rosh_hashana_greg_month = 9
    elif next_rosh_hashana_greg_day > 30 and pesach_greg_month == 4:
        next_rosh_hashana_greg_day = next_rosh_hashana_greg_day - 30
        next_rosh_hashana_greg_month = 10
        
    #print(next_rosh_hashana_greg_year, next_rosh_hashana_greg_month, next_rosh_hashana_greg_day)
    ############################################################################################################
    
    return (next_rosh_hashana_greg_year,next_rosh_hashana_greg_month,next_rosh_hashana_greg_day)

    
# פונקצייה שמחשבת כמה ימים עברו מאז ראש השנה העברי ועד היום
# היא ספציפית למיקרופייתון אך יכולה לעבוד בפייתון רגיל עם שינויים מתאימים לקבלת חותמת זמן
# פונקצייה זו משתמשת בפונקציות אחרות שהוגדרו למעלה
def get_days_from_rosh_hashana(greg_year, greg_month, greg_day):
     
    current_year = greg_year
    current_month = greg_month
    current_day = greg_day
    
    # הגדרת חותמת זמן של היום הנוכחי
    current_timestamp = get_mktime_date_only(current_year, current_month, current_day)
    
    # חישוב התאריך הלועזי של ראש השנה והגדרת חותמת זמן שלו
    rosh_hashana_greg = get_geus_rosh_hashana_greg(current_year)
    rosh_hashana_year, rosh_hashana_month, rosh_hashana_day = rosh_hashana_greg
    rosh_hashana_timestamp = get_mktime_date_only(rosh_hashana_year, rosh_hashana_month, rosh_hashana_day)
    
    # אם ראש השנה גדול מהיום הנוכחי כלומר שהוא עוד לא היה סימן שאנחנו צריכים את ראש השנה הקודם ולכן החישוב הוא על השנה הקודמת
    if rosh_hashana_timestamp > current_timestamp:
        # חישוב התאריך הלועזי של ראש השנה והגדרת חותמת זמן שלו
        rosh_hashana_greg = get_geus_rosh_hashana_greg(current_year-1) # הקטנת שנה
        rosh_hashana_year, rosh_hashana_month, rosh_hashana_day = rosh_hashana_greg
        rosh_hashana_timestamp = get_mktime_date_only(rosh_hashana_year, rosh_hashana_month, rosh_hashana_day)

      
    # חישוב ראש השנה הבא אחרי ראש השנה המבוקש
    next_rosh_hashana_greg = get_geus_rosh_hashana_greg(rosh_hashana_year+1) # חישוב ראש השנה הבא לאחר ראש השנה המבוקש 
    next_rosh_hashana_year, next_rosh_hashana_month, next_rosh_hashana_day = next_rosh_hashana_greg
    next_rosh_hashana_timestamp = get_mktime_date_only(next_rosh_hashana_year, next_rosh_hashana_month, next_rosh_hashana_day)

    # חישוב אורך השנה בימים
    length_heb_year_in_seconds = next_rosh_hashana_timestamp - rosh_hashana_timestamp
    length_heb_year_in_days = length_heb_year_in_seconds // (24 * 60 * 60)
    
    # חישוב הפרש הימים בין ראש השנה לבין היום
    days_from_rosh_hashana_in_seconds = current_timestamp - rosh_hashana_timestamp
    days_from_rosh_hashana = days_from_rosh_hashana_in_seconds // (24 * 60 * 60)
 
    rosh_hashana_heb_year_int = rosh_hashana_year + 3761 # זה בכוונה כך ולא 3760 כי מדובר על ראש השנה שחל לפני תחילת השנה הלועזית   
    
    return days_from_rosh_hashana, length_heb_year_in_days, rosh_hashana_heb_year_int

# פונקצייה שמחזירה את התאריך העברי הנוכחי כסטרינג וגם את מספר השנה העברית כאינט בהתבסס על הפונקציות הקודמות
def get_heb_date_and_holiday_from_greg_date(greg_year, greg_month, greg_day):
    days_from_rosh_hashana, length_heb_year_in_days, heb_year_int = get_days_from_rosh_hashana(greg_year, greg_month, greg_day)
    rosh_hashana_day, rosh_hashana_month = 1,1
    heb_day_int, heb_month_int = move_heb_date(rosh_hashana_day, rosh_hashana_month, length_heb_year_in_days, days_from_rosh_hashana)
    
    # האם השנה מעוברת
    is_leap_year = length_heb_year_in_days in [383, 384, 385]

    # האם כסלו חסר כלומר שיש בו רק 29 ימים זה תלוי באורך השנה
    is_kislev_29 = length_heb_year_in_days in [353, 383]
    
    # חישוב שם החודש והיום בעברית
    heb_day_string = heb_month_day_names(heb_day_int)
    heb_month_string = heb_month_names(heb_month_int, is_leap_year)
    heb_year_string = _num_to_str(heb_year_int, thousands=True, withgershayim=False)   
    heb_date_string = f"{heb_day_string}' {heb_month_string} {heb_year_string}"
    
    tuple_heb_date = (heb_day_int, heb_month_int, heb_year_int)
    
    holiday_name = get_holiday_name(heb_day_int, heb_month_int, is_leap_year)
    
    lite_holiday_name = get_lite_holiday_name(heb_day_int, heb_month_int, is_leap_year, is_kislev_29)

    is_rosh_chodesh = heb_day_int in [1,30]
    
    return heb_date_string, tuple_heb_date, holiday_name, lite_holiday_name, is_rosh_chodesh
    
# מחזיר תאריך עברי של היום הנוכחי כולל אפשרות ליום בשבוע
def get_today_heb_date_string(heb_week_day = False):
    # הגדרת הזמן הנוכחי המקומי מחותמת זמן לזמן רגיל
    tm = time.localtime(time.time())
    year, month, day, rtc_week_day, hour, minute, second, micro_second = (tm[0], tm[1], tm[2], tm[6], tm[3], tm[4], tm[5], 0)
    if heb_week_day:
        normal_weekday = get_normal_weekday(rtc_week_day)
        hebrew_weekday = heb_weekday_names(normal_weekday)
    heb_date_string, _, _, _, _, = get_heb_date_and_holiday_from_greg_date(year, month, day)
    return f'{hebrew_weekday}, {heb_date_string}' if heb_week_day else heb_date_string
    
def get_if_greg_is_heb_holiday(greg_year, greg_month, greg_day):
    _, _, holiday_name, lite_holiday_name, is_rosh_chodesh = get_heb_date_and_holiday_from_greg_date(greg_year, greg_month, greg_day)
    return holiday_name

def get_is_today_heb_holiday():
    year, month, day, rtc_week_day, hour, minute, second, micro_second = time.localtime(time.time())
    return get_if_greg_is_heb_holiday(year, month, day)

print(get_today_heb_date_string(heb_week_day=True))

##########################################################################################################################################
                           # 4 halacha_watch_C
##########################################################################################################################################

# ========================================================
# License: Personal Use Only  
# This file (`main_shemesh_s3.py`) is licensed for **personal use only**.  
# You may **not** modify, edit, or alter this file in any way.  
# Commercial use is strictly prohibited.  
# If you wish to use this file for business or organizational purposes,  
# please contact the author.  
# ========================================================

# משתנה גלובלי שמציין את גרסת התוכנה למעקב אחרי עדכונים
VERSION = "3/8/2025-C"

######################################################################################################################

# סיכום קצר על התוצאות המעשיות של הכפתורים בקוד הזה
# לחיצה על שתי הכפתורים בו זמנית כאשר מדליקים את המכשיר: עדכון תוכנת המכשיר
# אם כפתור המיקומים לחוץ בשעת הכניסה לתוכנה: עדכון השעון מהרשת 

# הבסיס של החישובים הגיע מכאן
# https://github.com/peterhinch/micropython-samples/tree/d2929df1b4556e71fcfd7d83afd9cf3ffd98fdac/astronomy
# לגבי בעיות עם המסך שבס"ד נפתרו ראו כאן
# https://github.com/Xinyuan-LilyGO/T-Display-S3/issues/300
# ובורא עולם תכנן שהפין שבחרתי להלחים עבור השעון החיצוני שזה פין 17 נשאר פעיל גם במהלך שינה עמוקה. זה טוב מאוד כי השעון החיצוני צריך למדוד טמפרטורה בשביל להיות מדוייק

import time, math, os, platform
import gc # חשוב נורא לניקוי הזיכרון
#from sun_moon_sgb import RiSet  # ספריית חישובי שמש
#from moonphase_sgb import MoonPhase  # ספריית חישובי שלב הירח
#import mpy_heb_date # לחישוב תאריך עברי מתאריך לועזי. ספרייה שלי
import tkinter as tk
from tkinter import font
from tkinter import messagebox

# לצורך לחיצה מדומה על המקלדת כדי שהמסך לא ייכבה אוטומטית
from pynput.keyboard import Controller as keyboard_Controller 
from pynput.keyboard import Key as keyboard_Key

is_windows = platform.system() == "Windows"
####################################################################3

# משתנה למעקב אחר מצב הכוח כלומר האם המכשיר כבוי או פועל
# המשמעות של זה מגיעה לידי ביטוי בפונקצייה הראשית: main
# זה משתנה חשוב נורא
power_state = True


# פונקצייה לקבלת הפרש השעות המקומי מגריניץ בלי התחשבות בשעון קיץ
# יש אפשרות להגדיר טרו או פאלס האם זה שעון קיץ או לא. כברירת מחדל זה לא
# אפשר להגדיר האם רוצים את ההפרש בשעות או בשניות וברירת המחדל היא בשעות
def get_generic_utc_offset(longitude_degrees, dst=False, in_seconds = False):
    offset = abs(round(longitude_degrees/15)) % 24
    offset = -offset if longitude_degrees < 0 else offset
    offset = offset + 1 if dst else offset
    return offset * 3600 if in_seconds else offset
    
# פונקצייה להפיכת טקסט כדי שעברית תהיה משמאל לימין
def reverse(string):
    return string if is_windows else "".join(reversed(string))

# מקבל מספר יום בשבוע לפי הנורמלי ומחזיר את מספר היום בשבוע לפי ההגדרות ב RTC
def get_rtc_weekday(weekday):
    rtc_weekday_dict = {1:6,2:0,3:1,4:2,5:3,6:4,7:5}
    return rtc_weekday_dict.get(weekday)

# מוציא את מספר היום בשבוע הנורמלי לפי סדר מתוך שעון המכשיר שמוגדר RTC
def get_normal_weekday(rtc_weekday):
    weekday_dict = {6:1,0:2,1:3,2:4,3:5,4:6,5:7}
    return weekday_dict.get(rtc_weekday)

# פונקצייה שמחזירה נכון או לא נכון האם כרגע נוהג שעון קיץ בישראל
# היא מתבססת על מה השעה והתאריך ברגע זה בשעון הפנימי של המיקרו בקר ולכן חייבים להגדיר אותו לפני שקוראים לפונקצייה זו
# שעון הקיץ מופעל בישראל בין יום שישי שלפני יום ראשון האחרון של חודש מרץ בשעה 02:00, לבין יום ראשון האחרון של חודש אוקטובר בשעה 02:00.
# השעה 2 בלילה של שינוי השעון הם בשעון ישראל ואילו שעון הבקר מוגדר לאיזור זמן גריניץ לכן הקדמתי בשעתיים לשעה 0 שעון גריניץ
def is_now_israel_DST():
    # קבלת השנה הנוכחית
    current_year = time.localtime()[0]
    
    # חישוב יום ראשון האחרון של מרץ
    march_last_sunday = time.mktime((current_year, 3, 31, 0, 0, 0, 0, 0, 0))
    while time.localtime(march_last_sunday)[6] != get_rtc_weekday(1):
        march_last_sunday -= 86400  # מורידים יום
    
    # חישוב יום שישי שלפני יום ראשון האחרון של מרץ
    # אם יום ראשון האחרון הוא ה-31, אז יום שישי לפניו יהיה ה-29.
    last_friday_march = march_last_sunday - 2 * 86400  # מורידים 2 ימים (שישי)

    # חישוב יום ראשון האחרון של אוקטובר
    october_last_sunday = time.mktime((current_year, 10, 31, 0, 0, 0, 0, 0, 0))
    while time.localtime(october_last_sunday)[6] != get_rtc_weekday(1): 
        october_last_sunday -= 86400  # מורידים יום
    
    # השוואה בין הזמן הנוכחי לתאריכים של שעון קיץ
    current_time = time.mktime(time.localtime())
    
    # שעון קיץ פעיל בין יום שישי שלפני יום ראשון האחרון של מרץ ועד יום ראשון האחרון של אוקטובר
    if last_friday_march <= current_time < october_last_sunday:
        return True  # שעון קיץ פעיל
    else:
        return False  # לא פעיל


# פונקצייה להמרת זמן מ-שניות ל- סטרינג שעות דקות ושניות, או רק ל- סטרינג דקות ושניות שבניתי בסיוע רובי הבוט
def convert_seconds(seconds, to_hours=False):        
    # חישוב מספר הדקות והשניות שיש בשעה אחת, והדפסתם בפורמט של דקות ושניות
    if to_hours:
        return f'{seconds // 3600 :02.0f}:{(seconds % 3600) // 60 :02.0f}:{seconds % 60 :02.0f}'
    else:
        return f'{seconds // 60 :02.0f}:{seconds % 60 :02.0f}'


# פונקצייה לפריטת שעה בשבר עשרוני לשניות
# אחר כך אפשר להשתמש בפונקצייה הקודמת להמרת של השניות לשעות דקות ושניות בפורמט שעון
# פונקצייה זו יכולה להיות שימושים לצורך נתינת עלייה ישרה בפורמט של שעות דקות ושניות אך כרגע אין לה שימוש בקוד זה
def decimal_hours_to_seconds(decimal_hours):
    hours = int(decimal_hours)  # קבלת השעות השלמות
    minutes_decimal = (decimal_hours - hours) * 60  # המרת החלק השברי לדקות
    minutes = int(minutes_decimal)  # קבלת הדקות השלמות
    seconds = (minutes_decimal - minutes) * 60  # המרת החלק השברי של הדקות לשניות
    total_seconds = hours * 3600 + minutes * 60 + seconds
    return total_seconds


#################################################################################################


# פונקצייה שמחשבת מה השעה הזמנית הנוכחית בהינתן הזמן הנוכחי וזמן הזריחה והשקיעה הקובעים
# כל הזמנים צריכים להינתן בפורמט חותמת זמן
# פונקצייה זו יכולה לפעול גם בכל פייתון רגיל היא לגמרי חישובית ולא תלוייה בכלום חוץ מהמשתנים שלה
def calculate_temporal_time(timestamp, sunrise_timestamp, sunset_timestamp):
    
        # בדיקה האם זה יום לפי בדיקה האם הזמן הנוכחי גדול משעת הזריחה וקטן משעת השקיעה
        is_day = timestamp >= sunrise_timestamp and timestamp < sunset_timestamp
        
        # חישוב מספר השקיעה מהזריחה לשקיעה אם זה יום או מהשקיעה לזריחה אם זה לילה
        day_or_night_length_seconds = sunset_timestamp - sunrise_timestamp if is_day else sunrise_timestamp - sunset_timestamp
        
        # חישוב מספר השניות בשעה זמנית אחת של היום לפי חלוקת אורך היום או הלילה ל 12
        seconds_per_temporal_hour_in_day_or_night = day_or_night_length_seconds / 12
        
        # חישוב כמה שניות עברו מאז הזריחה או השקיעה עד הזמן הנוכחי 
        time_since_last_sunrise_or_sunset = timestamp - (sunrise_timestamp if is_day else sunset_timestamp)
        
        # המרת השניות לפורמט שעות, דקות ושניות
        A = (time_since_last_sunrise_or_sunset / seconds_per_temporal_hour_in_day_or_night) + 0.0000001
        zmanit_hour = int(A)
        B = ((A - zmanit_hour) * 60) + 0.0000001
        zmanit_minute = int(B)
        C = ((B - zmanit_minute) * 60) + 0.0000001
        zmanit_second = int(C)

        # הדפסת השעה הזמנית המתאימה בפורמט שעות:דקות:שניות
        temporal_time = f'{zmanit_hour:02.0f}:{zmanit_minute:02.0f}:{zmanit_second:02.0f}'
        
        return temporal_time, seconds_per_temporal_hour_in_day_or_night


##############################################################################################################        
# הגדרת שמות עבור משתנים גלובליים ששומרים את כל הזמנים הדרושים לחישובים
sunrise, sunset, mga_sunrise, mga_sunset, yesterday_sunset, mga_yesterday_sunset, tomorrow_sunrise, mga_tomorrow_sunrise = [None] * 8
##############################################################################################################    

def get_sunrise_sunset_timestamps(current_timestamp, is_gra = True):
    
    # הצהרה על משתנים גלובליים ששומרים את הזמנים הדרושים
    global sunrise, sunset, mga_sunrise, mga_sunset, yesterday_sunset, mga_yesterday_sunset, tomorrow_sunrise, mga_tomorrow_sunrise
    
    #  חותמת זמן של רגע הזריחה והשקיעה היום
    # חוממת זמן של תחילת וסוף הדמדומים של מגן אברהם כרגע מוגדר לעיל מינוס 16
    sunrise_timestamp = sunrise if is_gra else mga_sunrise
    sunset_timestamp = sunset if is_gra else mga_sunset
           
    # בדיקה האם זה יום לפי בדיקה האם הזמן הנוכחי גדול משעת הזריחה וקטן משעת השקיעה
    is_day = current_timestamp >= sunrise_timestamp and current_timestamp < sunset_timestamp 
    
    #print("is_gra", is_gra)
    #print("is_day", is_day)
    
    if is_day:
        
        return sunrise_timestamp, sunset_timestamp
                
    else:
        # אם מדובר אחרי 12 בלילה ולפני הזריחה
        if current_timestamp < sunrise_timestamp:
            
            # הגדרת הזמן על אתמול וחישוב השקיעה של אתמול
            yesterday_sunset_timestamp = yesterday_sunset if is_gra else mga_yesterday_sunset
        
            return sunrise_timestamp, yesterday_sunset_timestamp
            
            
        # אם מדובר אחרי השקיעה ולפני השעה 12 בלילה
        elif (current_timestamp > sunrise_timestamp) and (current_timestamp >= sunset_timestamp):
            
            # הגדרת הזמן על מחר וחישוב הזריחה של מחר
            tomorrow_sunrise_timestamp = tomorrow_sunrise if is_gra else mga_tomorrow_sunrise
            
            return tomorrow_sunrise_timestamp, sunset_timestamp

################################################################################################################3

# כל המיקומים. כל מיקום מוגדר כמילון
# המיקום הראשון ברשימה יהיה ברירת המחדל אם לא מצליחים לקרוא מהקובץ שקובע מהו מיקום ברירת המחדל


locations = [
    
    {'heb_name': 'משווה-0-0', 'lat': 0.0, 'long': 0.0, 'altitude': 0.0, 'utc_offset': '', 'name': 'Equals 0-0'} , # קו המשווה בכוונה נמצא פעמיים כדי שבהעברת מיקומים תמיד יראו אותו ראשון
    {'heb_name': 'קו-המשווה', 'lat': 0.0, 'long': 0.0, 'altitude': 0.0, 'utc_offset': '', 'name': 'Equals 0-0'} ,
    {'heb_name': 'הקוטב-הצפוני', 'lat': 90.0, 'long': 0.0, 'altitude': 0.0, 'utc_offset': '', 'name': 'North Pole'} ,
    {'heb_name': 'ניו-יורק-ארהב', 'lat': 40.7143528, 'long': -74.0059731, 'altitude': 9.775694, 'utc_offset': '', 'name': 'New York US'} ,
    {'heb_name': 'אופקים', 'lat': 31.309, 'long': 34.61, 'altitude': 170.0, 'utc_offset': '', 'name': 'Ofakim IL'} ,
    {'heb_name': 'אילת', 'lat': 29.55, 'long': 34.95, 'altitude': 0.0, 'utc_offset': '', 'name': 'Eilat IL'} ,
    {'heb_name': 'אלעד', 'lat': 32.05, 'long': 34.95, 'altitude': 150.0, 'utc_offset': '', 'name': 'Elad IL'} ,
    {'heb_name': 'אריאל', 'lat': 32.10, 'long': 35.17, 'altitude': 521.0, 'utc_offset': '', 'name': 'Ariel IL'} ,
    {'heb_name': 'אשדוד', 'lat': 31.79, 'long': 34.641, 'altitude': 0.0, 'utc_offset': '', 'name': 'Ashdod IL'} ,
    {'heb_name': 'אשקלון', 'lat': 31.65, 'long': 34.56, 'altitude': 60.0, 'utc_offset': '', 'name': 'Ashkelon IL'} ,
    {'heb_name': 'באר-שבע', 'lat': 31.24, 'long': 34.79, 'altitude': 0.0, 'utc_offset': '', 'name': 'Beer Sheva IL'} ,
    {'heb_name': 'בית-שאן', 'lat': 32.5, 'long': 35.5, 'altitude': -120.0, 'utc_offset': '', 'name': 'Beit Shean IL'} ,
    {'heb_name': 'בית-שמש', 'lat': 31.74, 'long': 34.98, 'altitude': 300.0, 'utc_offset': '', 'name': 'Beit Shemesh IL'} ,
    {'heb_name': 'ביתר-עילית', 'lat': 31.69, 'long': 35.12, 'altitude': 800.0, 'utc_offset': '', 'name': 'Beitar Illit IL'} ,
    {'heb_name': 'בני-ברק', 'lat': 32.083156, 'long': 34.832722, 'altitude': 0.0, 'utc_offset': '', 'name': 'Bnei Brak IL'} ,
    {'heb_name': 'דימונה', 'lat': 31.07, 'long': 35.03, 'altitude': 570.0, 'utc_offset': '', 'name': 'Dimona IL'} ,
    {'heb_name': 'הר רומם', 'lat': 30.5100176, 'long': 34.6089109, 'altitude': 1000.0, 'utc_offset': '', 'name': 'Mount Romem IL'} ,
    {'heb_name': 'הרצליה', 'lat': 32.16, 'long': 34.84, 'altitude': 0.0, 'utc_offset': '', 'name': 'Herzliya IL'} ,
    {'heb_name': 'זכרון יעקב', 'lat': 32.57, 'long': 34.95, 'altitude': 170.0, 'utc_offset': '', 'name': 'Zichron Yaakov IL'} ,
    {'heb_name': 'חברון', 'lat': 31.53, 'long': 35.09, 'altitude': 950.0, 'utc_offset': '', 'name': 'Hebron IL'} ,
    {'heb_name': 'חדרה', 'lat': 32.43, 'long': 34.92, 'altitude': 53.0, 'utc_offset': '', 'name': 'Hadera IL'} ,
    {'heb_name': 'חיפה', 'lat': 32.8, 'long': 34.991, 'altitude': 300.0, 'utc_offset': '', 'name': 'Haifa IL'} ,
    {'heb_name': 'חרשה', 'lat': 31.944738, 'long': 35.1485598, 'altitude': 760.0, 'utc_offset': '', 'name': 'Harasha'} ,
    {'heb_name': 'טבריה', 'lat': 32.79, 'long': 35.531, 'altitude': 0.0, 'utc_offset': '', 'name': 'Tiberias IL'} ,
    {'heb_name': 'טלזסטון', 'lat': 31.78, 'long': 35.1, 'altitude': 720.0, 'utc_offset': '', 'name': 'Telzstone IL'} ,
    {'heb_name': 'ירוחם', 'lat': 30.99, 'long': 34.91, 'altitude': 0.0, 'utc_offset': '', 'name': 'Yeruham IL'} ,
    {'heb_name': 'ירושלים', 'lat': 31.776812, 'long': 35.235694, 'altitude': 750.0, 'utc_offset': '', 'name': 'Jerusalem IL'} ,
    {'heb_name': 'כוכב השחר', 'lat': 31.96, 'long': 35.34, 'altitude': 577.0, 'utc_offset': '', 'name': 'Cochav Hashachar IL'} ,
    {'heb_name': 'כרמיאל', 'lat': 32.915, 'long': 35.292, 'altitude': 315.0, 'utc_offset': '', 'name': 'Carmiel IL'} ,
    {'heb_name': 'לוד', 'lat': 31.95, 'long': 34.89, 'altitude': 0.0, 'utc_offset': '', 'name': 'Lod IL'} ,
    {'heb_name': 'מגדל-העמק', 'lat': 32.67, 'long': 35.23, 'altitude': 0.0, 'utc_offset': '', 'name': 'Migdal Haemek IL'} ,
    {'heb_name': 'מודיעין-עילית', 'lat': 31.940826, 'long': 35.037057, 'altitude': 320.0, 'utc_offset': '', 'name': "Modi'in Illit IL"} ,
    {'heb_name': 'מיצד', 'lat': 31.585503, 'long': 35.187587, 'altitude': 937.0, 'utc_offset': '', 'name': 'Meizad IL'} ,
    {'heb_name': 'מירון', 'lat': 32.98, 'long': 35.43, 'altitude': 700.0, 'utc_offset': '', 'name': 'Miron IL'} ,
    {'heb_name': 'מצפה רמון', 'lat': 30.6097894, 'long': 34.8120107, 'altitude': 855.0, 'utc_offset': '', 'name': 'Mitzpe Ramon IL'} ,
    {'heb_name': 'נהריה', 'lat': 33.01, 'long': 35.1, 'altitude': 25.0, 'utc_offset': '', 'name': 'Nahariya IL'} ,
    {'heb_name': 'נחליאל', 'lat': 31.9743, 'long': 35.14038, 'altitude': 575.0, 'utc_offset': '', 'name': 'Nahaliel'} ,
    {'heb_name': 'נצרת-עילית', 'lat': 32.7, 'long': 35.32, 'altitude': 0.0, 'utc_offset': '', 'name': 'Nazareth Illit IL'} ,
    {'heb_name': 'נתיבות', 'lat': 31.42, 'long': 34.59, 'altitude': 142.0, 'utc_offset': '', 'name': 'Netivot IL'} ,
    {'heb_name': 'נתניה', 'lat': 32.34, 'long': 34.86, 'altitude': 0.0, 'utc_offset': '', 'name': 'Netanya IL'} ,
    {'heb_name': 'עכו', 'lat': 32.93, 'long': 35.08, 'altitude': 0.0, 'utc_offset': '', 'name': 'Ako IL'} ,
    {'heb_name': 'עמנואל', 'lat': 32.16, 'long': 35.13, 'altitude': 406.0, 'utc_offset': '', 'name': 'Emmanuel IL'} ,
    {'heb_name': 'עפולה', 'lat': 32.6, 'long': 35.29, 'altitude': 60.0, 'utc_offset': '', 'name': 'Afula IL'} ,
    {'heb_name': 'ערד', 'lat': 31.26, 'long': 35.21, 'altitude': 640.0, 'utc_offset': '', 'name': 'Arad IL'} ,
    {'heb_name': 'פתח-תקווה', 'lat': 32.09, 'long': 34.88, 'altitude': 0.0, 'utc_offset': '', 'name': 'Petah Tikva IL'} ,
    {'heb_name': 'צפת', 'lat': 32.962, 'long': 35.496, 'altitude': 850.0, 'utc_offset': '', 'name': 'Zefat IL'} ,
    {'heb_name': 'קצרין', 'lat': 32.98, 'long': 35.69, 'altitude': 0.0, 'utc_offset': '', 'name': 'Katzrin IL'} ,
    {'heb_name': 'קרית-גת', 'lat': 31.61, 'long': 34.77, 'altitude': 159.0, 'utc_offset': '', 'name': 'Kiryat Gat IL'} ,
    {'heb_name': 'קרית-שמונה', 'lat': 33.2, 'long': 35.56, 'altitude': 0.0, 'utc_offset': '', 'name': 'Kiryat Shmona IL'} ,
    {'heb_name': 'ראש-העין', 'lat': 32.08, 'long': 34.95, 'altitude': 90.0, 'utc_offset': '', 'name': 'Rosh HaAyin IL'} ,
    {'heb_name': 'ראשון-לציון', 'lat': 31.96, 'long': 34.8, 'altitude': 0.0, 'utc_offset': '', 'name': 'Rishon Lezion IL'} ,
    {'heb_name': 'רחובות', 'lat': 31.89, 'long': 34.81, 'altitude': 76.0, 'utc_offset': '', 'name': 'Rechovot IL'} ,
    {'heb_name': 'רכסים', 'lat': 32.74, 'long': 35.08, 'altitude': 154.0, 'utc_offset': '', 'name': 'Rechasim IL'} ,
    {'heb_name': 'רמלה', 'lat': 31.92, 'long': 34.86, 'altitude': 0.0, 'utc_offset': '', 'name': 'Ramla IL'} ,  
    {'heb_name': 'רעננה', 'lat': 32.16, 'long': 34.85, 'altitude': 71.0, 'utc_offset': '', 'name': 'Raanana IL'} ,
    {'heb_name': 'שדרות', 'lat': 31.52, 'long': 34.59, 'altitude': 0.0, 'utc_offset': '', 'name': 'Sderot IL'} ,
    {'heb_name': 'שילה', 'lat': 32.05, 'long': 35.29, 'altitude': 719.0, 'utc_offset': '', 'name': 'Shilo IL'} ,
    {'heb_name': 'תל-אביב-חולון', 'lat': 32.01, 'long': 34.75, 'altitude': 0.0, 'utc_offset': '', 'name': 'Tel Aviv-Holon IL'} ,
    {'heb_name': 'תפרח', 'lat': 31.32, 'long': 34.67, 'altitude': 173.0, 'utc_offset': '', 'name': 'Tifrach IL'} ,
    
    {'heb_name': 'אומן-אוקראינה', 'lat': 48.74732, 'long': 30.23332, 'altitude': 211.0, 'utc_offset': '', 'name': 'Uman UA'} ,
    {'heb_name': 'אמסטרדם-הולנד', 'lat': 52.38108, 'long': 4.88845, 'altitude': 15.0, 'utc_offset': '', 'name': 'Amsterdam NL'} ,
    {'heb_name': 'וילנא-ליטא', 'lat': 54.672298, 'long': 25.2697, 'altitude': 112.0, 'utc_offset': '', 'name': 'Vilnius LT'} ,
    {'heb_name': "ז'שוב-פולין", 'lat': 50.0332, 'long': 21.985848, 'altitude': 209.0, 'utc_offset': '', 'name': 'Rzeszow PL'} ,
    {'heb_name': 'טרולהטן-שבדיה', 'lat': 58.28, 'long': 12.28, 'altitude': 28.0, 'utc_offset': '', 'name': 'Trollhattan SE'} ,  
    {'heb_name': 'לונדון-אנגליה', 'lat': 51.5001524, 'long': -0.1262362, 'altitude': 14.605533, 'utc_offset': '', 'name': 'London GB'} ,
    {'heb_name': 'לייקווד-ארהב', 'lat': 40.07611, 'long': -74.21993, 'altitude': 16.0, 'utc_offset': '', 'name': 'Lakewood US'} ,
    {'heb_name': 'מוסקווה-רוסיה', 'lat': 55.755786, 'long': 37.617633, 'altitude': 151.189835, 'utc_offset': '', 'name': 'Moscow RU'} ,
    {'heb_name': 'סטוקהולם-שבדיה', 'lat': 59.33, 'long': 18.06, 'altitude': 28.0, 'utc_offset': '', 'name': 'Stockholm SE'} ,
    {'heb_name': "פראג-צ'כיה", 'lat': 50.0878114, 'long': 14.4204598, 'altitude': 191.103485, 'utc_offset': '', 'name': 'Prague CZ'} ,
    {'heb_name': 'פריז-צרפת', 'lat': 48.8566667, 'long': 2.3509871, 'altitude': 0.0, 'utc_offset': '', 'name': 'Paris FR'} ,
    {'heb_name': 'פרנקפורט-גרמניה', 'lat': 50.1115118, 'long': 8.6805059, 'altitude': 106.258285, 'utc_offset': '', 'name': 'Frankfurt DE'} ,
    {'heb_name': 'קהיר-מצרים', 'lat': 30.00022, 'long': 31.231873, 'altitude': 23.0, 'utc_offset': '', 'name': 'Cairo EG'} ,
    {'heb_name': 'רומא-איטליה', 'lat': 41.8954656, 'long': 12.4823243, 'altitude': 19.704413, 'utc_offset': '', 'name': 'Rome IT'} ,
    {'heb_name': 'רמרוג-צרפת רת', 'lat': 48.518606, 'long': 4.3034152, 'altitude': 101.0, 'utc_offset': '', 'name': 'Ramerupt FR'} ,
    
       
    ]



#  ההסברים מורכבים משני חלקים כל אחד: הסבר וערך. ההסבר עובר בסוף רוורס ולכן אם יש בו מספרים חייבים לעשות להם רוורס כאן כדי שהרוורס הסופי יישר אותם 
esberim = [
    
        ["ליציאה לחצו במקלדת על",f"Esc"],
        ["או לחצו בעכבר על לחצן אמצעי",""],
        ["מעבר בין מיקומים בחיצים ימין ושמאל", ""],
        ["או באמצעות גלילה בעכבר", ""],
        ["שיפט+חץ ימין: קביעת מיקום ברירת מחדל", ""],
        ["חץ למעלה: חזרה למיקום ברירת מחדל", ""],
        
        ["שעון ההלכה גרסה",f"{VERSION}"],
        [" מאת: שמחה גרשון בורר - כוכבים וזמנים",""],
        [reverse("sgbmzm@gmail.com"), ""],
        ["כל הזכויות שמורות - להלן הסברים", ""],
        
        [f"כשהשעון מכוון: דיוק הזמנים {reverse('10')} שניות", ""],
        
        [" התאריך העברי מתחלף בשקיעה", ""],
        
        [" מתחת גרא/מגא:  דקות בשעה זמנית", ""],
        [" מתחת שמש/ירח:  אזימוט שמש/ירח", ""],
        ["אזימוט = מעלות מהצפון, וכדלהלן", ""],
        [f"צפון={reverse('0/360')}, מז={reverse('90')}, ד={reverse('180')}, מע={reverse('270')}", ""],
        ["  שלב הירח במסלולו החודשי - באחוזים", ""],
        [f"מולד={reverse('0/100')}, ניגוד={reverse('50')}, רבעים={reverse('25/75')}", ""],
    
        ["רשימת זמני היום בשעות זמניות", ""],
        ["זריחה ושקיעה במישור", "00:00"],
        ["סוף שמע ביום/רבע הלילה", "03:00"],
        ["  סוף תפילה ביום/שליש הלילה", "04:00"],
        ["חצות יום ולילה", "06:00"],
        ["מנחה גדולה", "06:30"],
        ["מנחה קטנה", "09:30"],
        ["פלג המנחה", "10:45"],
        [f"מגא מחושב לפי {reverse('-16°')} בבוקר ובערב", ""],
        
        ["   זמנים במעלות כשהשמש תחת האופק", ""],
        ["זריחת ושקיעת מרכז השמש", "0.0°"],
        ["  זריחה ושקיעה במישור", "-0.833°"],
        
        [f"זמני צאת הכוכבים {reverse('3/4')} מיל במעלות", ""],
        [f"לפי מיל של {reverse('18')} דקות", "-3.65°"],
        [f"לפי מיל של {reverse('22.5')} דקות", "-4.2°"],
        [f"לפי מיל של {reverse('24')} דקות", "-4.61°"],
        ["צאת כוכבים קטנים רצופים", "-6.3°"],
        
        ["  מעלות: עלות השחר/צאת כוכבים דרת", ""],
        [f"לפי 4 מיל של {reverse('18')} דקות", "-16.02°"],
        [f"לפי 4 מיל של {reverse('22.5')} דקות", "-19.75°"],
        [f"לפי 5 מיל של {reverse('24')} דקות", "-25.8°"],
        ["משיכיר/תחילת ציצית ותפילין", "-10.5°"],
        
        
        ["זמנים נוספים", ""],
        ["להימנע מסעודה בערב שבת", "09:00"],
        ["סוף אכילת חמץ", "04:00"],
        ["סוף שריפת חמץ", "05:00"],

        
        ["להלן תנאי מינימום לראיית ירח ראשונה", ""],
        [f"שלב {reverse('3%')}; והפרש גובה שמש-ירח {reverse('8°')}", ""],
    
    ]

# פונקצייה שמחזירה את השעה במיקום שבו נמצאים כרגע כחותמת זמן
def get_current_location_timestamp(manual_time = False):
    
    # הגדרות מאוד חשובות על איזה זמן יתבצעו החישובים
    # בתחילת הקוד גרמנו שהשעון החיצוני וגם הפנימי מעודכנים בשעה בגריניץ כלומר באיזור זמן UTC-0 . כעת צריך להמיר לשעון מקומי במיקום הנוכחי
    rtc_system_timestamp =  time.time() # או: time.mktime(time.localtime())
    #rtc_system_timestamp = time.mktime((2025, 7, 25, 20, 5, 45, 0, 0,-1)) # זה לבדיקה בלבד כשרוצים להזין זמן ידני
    current_utc_timestamp =  rtc_system_timestamp # כי בתחילת הקוד גרמנו שהשעון החיצוני יעדכן את השעון הפנימי בשעה באיזור זמן UTC-0
    # בדיקה האם המיקום הנוכחי הוא משווה 00 או הקוטב הצפוני אפס כי שם אני לא רוצה שיהיה שעון קיץ
    is_location_mashve_or_kotev = location["long"] == 0.0 and location["lat"] == 0.0 or location["long"] == 0.0 and location["lat"] == 90.0
    is_location_dst = True if is_now_israel_DST() and not is_location_mashve_or_kotev else False # כרגע כל שעון הקיץ או לא שעון קיץ נקבע לפי החוק בישראל גם עבור מקומות אחרים
    location_offset_hours = get_generic_utc_offset(location["long"], dst=is_location_dst) # חישוב הפרש הזמן מגריניץ עבור המיקום הנוכחי בשעות
    location_offset_seconds = get_generic_utc_offset(location["long"], dst=is_location_dst, in_seconds = True) # חישוב הפרש הזמן בשניות
    current_location_timestamp = current_utc_timestamp + location_offset_seconds # חותמת הזמן המקומית היא UTC-0 בתוספת הפרש השניות המקומי
    # עכשיו הגענו לנתון הכי חשוב שהוא חותמת הזמן המקומית הנוכחית
    return current_location_timestamp, location_offset_hours, location_offset_seconds



########################################################################################3
# פונקציות לטיפול במיקומים
# פונקצייה שמחזירה את המיקום במחשב שבו שומרים את הקובץ של מיקום ברירת המחדל
def get_settings_path():
    
    APP_NAME = "Halacha_watch"
    
    if platform.system() == "Windows":
        base_dir = os.getenv("LOCALAPPDATA")
    else:
        base_dir = os.path.expanduser("~/.config")

    config_dir = os.path.join(base_dir, APP_NAME)
    os.makedirs(config_dir, exist_ok=True)
    return os.path.join(config_dir, "hw_default_location_index.txt")

# פונקצייה לשמירת מיקום ברירת המחדל
def save_default_location_index(event=None):
    # הכרזה על משתנים גלובליים שיטופלו בלחיצה על הכפתור
    global location_index
    global location
    path = get_settings_path()
    answer = messagebox.askyesno(reverse("אישור"), reverse(f"האם לשמור את '{location['heb_name']}' כמיקום ברירת מחדל?"))
    if answer:
        try:
            with open(path, "w") as f:
                f.write(str(location_index))
                messagebox.showinfo(reverse("אישור"), reverse(f"מיקום ברירת מחדל הוגדר בהצלחה: {location['heb_name']}"))
        except Exception as e:
            messagebox.showinfo(reverse("שגיאה"), f"{e}")

# פונקצייה לקריאת מיקום ברירת המחדל מתוך הקובץ
def load_default_location_index():
    path = get_settings_path()
    if os.path.exists(path):
        with open(path, "r") as f:
            return int(f.read())
    else:
        return 0  # ברירת מחדל אם הקובץ לא קיים

# פונקצייה להחלפת מיקום
def switch_location(event=None):
    global location, location_index
    # event.keysym מיועד לאירועי מקלדת event.num מיועד לאירועי עכבר ו event.delta מיועד לגלילת עכבר בווינדוס
    if event.keysym == "Right" or event.num in [3,4] or event.delta > 0: 
        location_index = (location_index + 1) % len(locations)
    elif event.keysym == "left" or event.num in [1,5] or event.delta < 0:
        location_index = (location_index - 1) % len(locations)

    location = locations[location_index]
    

# פונקצייה שמקפיצה בחזרה את התוכנה להיות על מיקום ברירת המחדל
def go_to_default_location(event=None):
    # הצהרה על משתנים גלובליים
    global location, location_index
    # מחזיר את המיקום הנוכחי להיות מיקום ברירת מחדל
    default_index = load_default_location_index()
    location = locations[default_index] if 0 <= default_index < len(locations) else locations[0]
    # מאפס את המיקום שאוחזים בו בדפדוף ברשימת המיקומים כך שהדפדוף הבא יתחיל מהתחלה ולא מהמיקום האינדקסי של מיקום ברירת המחדל
    location_index = 0

    
# הגדרת משתנה גלובלי חשוב מאוד שקובע מה המיקום הנוכחי שעליו מתבצעים החישובים
# משתנה זה נקבע לפי המיקום האינדקסי ששמור בקובץ מיקום ברירת מחדל תוך בדיקה שהאינדקס לא חורג מגבולות הרשימה ואם כן חורג אז יוגדר המיקום האפס כברירת מחדל
# קריאת המיקום מתוך הרשימה בהתאם למספר שבקובץ
default_index = load_default_location_index()
location = locations[default_index] if 0 <= default_index < len(locations) else locations[0] 
# אינדקס המיקום הנוכחי משתנה גלובלי חשוב מאוד לצורך התקדמות ברשימת המיקומים
# הגדרה שלו על אפס גורמת שכל דפדוף ברשימת המיקומים יתחיל מהתחלה ולא ימשיך מהמיקום האינדקסי של מיקום ברירת המחדל
location_index = 0

##############################################################################################


# משתנה לשליטה על איזה נתונים יוצגו בהסברים במסך של שעון ההלכה בכל שנייה
current_screen_halach_clock = 0.0  #

# משתנה לשליטה אלו נתונים יוצגו בשורת הזמנים 
current_screen_zmanim = 0

# מונה לדעת מתי ללחוץ לחיצה ווירטואלית על שיפט כדי למנוע כיבוי מסך
counter_shift = 0.0

###########################################################################################################3

# הקמת החלון הגרפי

root_hw = tk.Tk()
root_hw.attributes('-fullscreen', True) # מסך מלא
root_hw.configure(bg='black') # רקע שחור 

root_hw.bind("<Escape>", lambda e: root_hw.destroy()) # כיבוי בלחיצה על אסקייפ
root_hw.bind("<Button-2>", lambda e: root_hw.destroy()) # כיבוי בלחיצה על לחצן אמצעי בעכבר

root_hw.bind("<Shift-Right>", save_default_location_index) # שיפט וחץ ימינה שומרים מיקום ברירת מחדל

root_hw.bind("<Right>", switch_location) # חץ ימינה מחליף מיקום
root_hw.bind("<Left>", switch_location) # חץ שמאלה מחליף מיקום
root_hw.bind("<Up>", go_to_default_location) # חץ למעלה מחזיר למיקום ברירת המחדל

#root_hw.bind("<Button-1>", switch_location) # לחצן שמאלי בעכבר מחליף מיקום
#root_hw.bind("<Button-3>", switch_location) # לחצן ימני בעכבר מחליף מיקום

# תמיכה גם בלינוקס וגם בווינדוס בגלילת העכבר
root_hw.bind("<MouseWheel>", switch_location)      # Windows scroll up & down
root_hw.bind("<Button-4>", switch_location)        # Linux scroll up
root_hw.bind("<Button-5>", switch_location)        # Linux scroll down

# עושה שלא יראו את העכבר על המסך של שעון ההלכה
root_hw.config(cursor="none")

'''
# כבר אין צורך בזה כי עברתי לשיטה טובה יותר וזה נשאר רק לזיכרון
if is_windows:
    ######################################################################################
    import ctypes
    # מונע מצב שינה וכיבוי מסך
    ES_CONTINUOUS = 0x80000000 # אומר: ההגדרה הזו תישאר בתוקף עד שנשנה אותה שוב או שהתהליך יסתיים.
    ES_DISPLAY_REQUIRED = 0x00000002 # אומר: שמור על המסך פעיל – אל תכבה אותו
    ctypes.windll.kernel32.SetThreadExecutionState(ES_CONTINUOUS | ES_DISPLAY_REQUIRED)
    ######################################################################################
'''

# הגדרות החלון
# רזולוציית המסך שבשימוש כרגע
screen_width = root_hw.winfo_screenwidth()
screen_height = root_hw.winfo_screenheight()

base_width = 320
base_height = 170
scale_x = screen_width / base_width
scale_y = screen_height / base_height
scale = min(scale_x, scale_y)

# פונקציה ליצירת גופנים עם קנה מידה
def scaled_font(name, size, weight="normal"):
    return (name, int(size * scale), weight)

canvas = tk.Canvas(root_hw, width=screen_width, height=screen_height, bg="black", highlightthickness=0)
canvas.pack()

# טקסטים ונתונים לשיבוץ בחלון הגרפי
hw_green = "lime"

# איזור כותרת
title_id = canvas.create_text(160 * scale, 12 * scale, text="", fill=hw_green, font=scaled_font("miriam", 12, "bold"))

# איזור תאריך עברי
heb_date_rect_id = canvas.create_rectangle(0, 21 * scale, screen_width, 38 * scale, fill="black")
heb_date_id = canvas.create_text(160 * scale, 30 * scale, text="", fill="white", font=scaled_font("miriam", 14, "bold"))
canvas.create_line(0, 45 * scale, screen_width, 45 * scale, fill="yellow")

# איזור שעה זמנית גרא ומגא
canvas.create_text(300 * scale, 55 * scale, text=reverse("גרא"), fill="white", font=scaled_font("miriam", 13))
minutes_in_gra_temporal_hour_id = canvas.create_text(300 * scale, 70 * scale, text="", fill="turquoise", font=scaled_font("miriam", 13))
gra_temporal_hour_id = canvas.create_text(210 * scale, 64 * scale, text="", fill=hw_green, font=scaled_font("miriam", 30, "bold"))

canvas.create_text(120 * scale, 55 * scale, text=reverse("מגא"), fill="white", font=scaled_font("miriam", 13))
minutes_in_mga_temporal_hour_id = canvas.create_text(120 * scale, 70 * scale, text="", fill="turquoise", font=scaled_font("miriam", 13))
mga_temporal_hour_id = canvas.create_text(53 * scale, 64 * scale, text="", fill=hw_green, font=scaled_font("miriam", 18, "bold"))
canvas.create_line(0, 80 * scale, screen_width, 80 * scale, fill="yellow")

# איזור מידע על שמש וירח
canvas.create_text(300 * scale, 90 * scale, text=reverse("שמש"), fill="white", font=scaled_font("miriam", 13))
sun_az_id = canvas.create_text(300 * scale, 108 * scale, text="", fill="turquoise", font=scaled_font("miriam", 13))
sun_alt_id = canvas.create_text(210 * scale, 102 * scale, text="", fill=hw_green, font=scaled_font("miriam", 30, "bold"))

canvas.create_text(120 * scale, 90 * scale, text=reverse("ירח"), fill="white", font=scaled_font("miriam", 13))
moon_az_id = canvas.create_text(120 * scale, 108 * scale, text="", fill="turquoise", font=scaled_font("miriam", 13))
moon_alt_id = canvas.create_text(53 * scale, 93 * scale, text="", fill=hw_green, font=scaled_font("miriam", 18, "bold"))
moon_phase_id = canvas.create_text(45 * scale, 111 * scale, text="", fill="turquoise", font=scaled_font("miriam", 15, "bold"))
canvas.create_line(0, 120 * scale, screen_width, 120 * scale, fill="yellow")

# איזור הסברים מתחלף
hesberim_id = canvas.create_text(160 * scale, 132 * scale, text="", fill="white", font=scaled_font("miriam", 14))
canvas.create_line(0, 143 * scale, screen_width, 143 * scale, fill="yellow")

# איזור שעון רגיל ותאריך לועזי ואיזור הזמן
utc_offset_id = canvas.create_text(270 * scale, 157 * scale, text="", fill="white", font=scaled_font("miriam", 18))
time_id = canvas.create_text(180 * scale, 157 * scale, text="", fill=hw_green, font=scaled_font("miriam", 20, "bold"))
greg_date_id = canvas.create_text(65 * scale, 157 * scale, text="", fill="white", font=scaled_font("miriam", 18))
canvas.create_line(0, 166 * scale, screen_width, 166 * scale, fill="yellow")

# איזור שקיים רק במחשב ולא בשעון ההלכה הפיזי שמציג זמני הלכה קשיחים בשעון רגיל
zmanim_id = canvas.create_text(160 * scale, 174 * scale, text="", fill="magenta", font=scaled_font("miriam", 10))

######################################################################################################################3

# הפונקצייה הראשית שבסוף גם מפעילה את הנתונים על המסך
def main_halach_clock():
    
     
    # קבלת הזמן המקומי למיקום המבוקש כחותמת זמן - באמצעות פונקצייה שהוגדרה לעיל    
    current_location_timestamp, location_offset_hours, location_offset_seconds = get_current_location_timestamp()
    current_timestamp = current_location_timestamp
            
    # משתנה ששולט על חישוב גובה השמש במעלות לשיטת המג"א ונועד במקור לחישוב דמדומים
    # אם כותבים 16 זה אומר מינוס 16
    # אם רוצים פלוס אז אולי צריך לעשות +16 אבל לא יודע אם זה יעבוד
    # אם עושים None או False או 0 זה לא מחושב כלל (ולכן אם במקרה רוצים כאן זריחה גיאומטרית חייבים להגדיר 0.00001)
    MGA_deg = 16 
    
    # הצהרה על משתנים גלובליים ששומרים את הזמנים הדרושים
    global sunrise, sunset, mga_sunrise, mga_sunset, yesterday_sunset, mga_yesterday_sunset, tomorrow_sunrise, mga_tomorrow_sunrise
    
    # ריקון כל המשתנים כדי שלא ישתמשו בנתונים לא נכונים
    sunrise, sunset, mga_sunrise, mga_sunset, yesterday_sunset, mga_yesterday_sunset, tomorrow_sunrise, mga_tomorrow_sunrise = [None] * 8
    
        
    # יצירת אובייקט RiSet
    RiSet.tim = round(current_timestamp) ############### אם לא מגדירים את זה אז הזמן הוא לפי הזמן הפנימי של הבקר
    #RiSet.sinho_sun_riset = 0.0 # אם רוצים שזריחה ושקיעה של השמש יהיו לפי זריחה ושקיעה גיאומטריים ולא לפי מינוס 0.833. אם לא מגדירים אז כברירת מחדל יהיה 0.833 מינוס
    riset = RiSet(lat=location["lat"], long=location["long"], lto=location_offset_hours, tl=MGA_deg) # lto=location_offset_hours
    
    ############# חישוב גובה ואזימוט של השמש והירח ברגע הנוכחי ###########
    ####### חובה לעשות את זה כאן כאשר ריסט מוגדר עדיין על היון הנוכחי כי בשלב הבא לפעמים מגדירים את ריסט על היום הבא או הקודם ######
    
    # הגדרת הזמן הנוכחי המקומי מחותמת זמן לזמן רגיל
    tm = time.gmtime(current_timestamp) # אסור להשתמש כאן ב time.localtime כי זה בפייתון רגיל מחזיר זמן מקומי של המחשב
    year, month, day, rtc_week_day, hour, minute, second, micro_second = (tm[0], tm[1], tm[2], tm[6], tm[3], tm[4], tm[5], 0)
        
    # חישוב מה השעה הנוכחית בשבר עשרוני
    current_hour = (hour + (minute / 60) + (second / 3600)) - location_offset_hours
    
    # חישוב גובה השמש והירח וגם אזימוט ועלייה ישרה (עלייה ישרה בשעות עשרוני) ברגע זה. כלומר בשעה הנוכחית בשבר עשרוני
    # לדעת את גובה השמש והירח אפשר גם במיקום שאין בו זריחות ושקיאות וזה לא מחזיר שגיאה אלא מחזיר None שזה כמו אפס
    s_alt, s_az, s_ra, s_dec = riset.alt_az_ra_dec(current_hour, sun=True)
    m_alt, m_az, m_ra, m_dec = riset.alt_az_ra_dec(current_hour, sun=False)
    
    ########## חישובי זריחות ושקיעות היום וגם אתמול או מחר הדרושים לחישוב שעון שעה זמנית #############
      
    # שמירת כל הנתונים על היום הנוכחי כי כולם נוצרים ביחד בעת הגדרת "riset" או בעת שמשנים לו יום
    sunrise, sunset, mga_sunrise, mga_sunset = riset.sunrise(1), riset.sunset(1), riset.tstart(1), riset.tend(1)
    
    # אם מדובר אחרי 12 בלילה ולפני הזריחה לפי אחת משתי השיטות ההלכתיות (ויודעים את זה לפי ששעת הזריחה מאוחרת מהרגע הנוכחי) השעה הזמנית מתחילה בשקיעה של אתמול
    # מגדרים את יום האתמול ושומרים את כל הנתונים הדרושים עכשיו או בעתיד על יום האתמול    
    # אם בעתיד ירצו שעון ערבי מהשקיעה הקודמת יהיו חייבים להפעיל את החישוב הזה בקביעות ולא רק אחרי 12 בלילה ולפני הזריחה
    if (sunrise and current_timestamp < sunrise) or (mga_sunrise and current_timestamp < mga_sunrise):
        riset.set_day(-1)
        yesterday_sunset, mga_yesterday_sunset = riset.sunset(1), riset.tend(1) if mga_sunrise else None
        tomorrow_sunrise, mga_tomorrow_sunrise = None, None # לא חייבים את זה אבל זה מוסיף לביטחות שלא יתבצעו חישובים על נתונים לא נכונים
        
    # אם מדובר אחרי השקיעה לפי אחת השיטות ולפני השעה 12 בלילה השעה הזמנית מסתיימת בזריחה של מחר
    # מגדירים את יום המחר ושומרים את כל הנתונים הדרושים עכשיו או בעתיד על יום המחר
    elif (sunrise and sunset and current_timestamp > sunrise and current_timestamp >= sunset) or (mga_sunrise and mga_sunset and current_timestamp > mga_sunrise and current_timestamp >= mga_sunset):
        riset.set_day(1)
        tomorrow_sunrise, mga_tomorrow_sunrise  = riset.sunrise(1), riset.tstart(1) if mga_sunrise else None 
        yesterday_sunset, mga_yesterday_sunset = None, None # לא חייבים את זה אבל זה מוסיף לביטחות שלא יתבצעו חישובים על נתונים לא נכונים
   
    # הדפסות לניסיון כשיש בעיות
    #print("mga_sunrise",time.strftime("%H:%M:%S %d/%m/%Y",time.gmtime(mga_sunrise)))
    #print("mga_sunset", time.strftime("%H:%M:%S %d/%m/%Y",time.gmtime(mga_sunset)))
    
       
    ################## חישוב השעה הזמנית הנוכחית גרא ומגא  ##################
   
    # כל החישובים נעשים רק אם יש זריחה ושקיעה ביממה זו במיקום זה והזריחה היא לפני השקיעה. כי אולי במיקום הזה אין בכלל זריחה ושקיעה ביום זה
    if sunrise and sunset and sunrise < sunset:
        
        # חישוב מה הם הזריחה והשקיעה הקובעים את השעון של שעה זמנית באמצעות פונקצייה שהוגדרה למעלה    
        sunrise_timestamp, sunset_timestamp = get_sunrise_sunset_timestamps(current_timestamp, is_gra = True)
         
        # חישוב שעון שעה זמנית על הזריחה והשקיעה באמצעות פונקצייה שהוגדרה למעלה
        temporal_time, seconds_in_temporal_hour = calculate_temporal_time(current_timestamp, sunrise_timestamp, sunset_timestamp)
        minutes_in_temporal_hour = str(round(seconds_in_temporal_hour / 60,1)) # str(convert_seconds(seconds_in_temporal_hour))
             
    else:
        
        temporal_time = reverse("שגיאה  ")
        minutes_in_temporal_hour = 0.0
        
       
    # כל החישובים נעשים רק אם יש זריחה ושקיעה של מגא ביממה זו במיקום זה והזריחה היא לפני השקיעה.
    # כי אולי במיקום הזה אין בכלל זריחה ושקיעה ביום זה כלומר שהשמש לא יורדת אל מתחת האופק בצורה מסודרת ביממה זו
    if mga_sunrise and mga_sunset and mga_sunrise < mga_sunset:

        # חישוב מחדש עבור שיטת מגן אברהם    
        # חישוב מה הם הזריחה והשקיעה הקובעים את השעון של שעה זמנית באמצעות פונקצייה שהוגדרה למעלה    
        mga_sunrise_timestamp, mga_sunset_timestamp = get_sunrise_sunset_timestamps(current_timestamp, is_gra = False)
         
        # חישוב שעון שעה זמנית על הזריחה והשקיעה באמצעות פונקצייה שהוגדרה למעלה
        mga_temporal_time, seconds_in_mga_temporal_hour = calculate_temporal_time(current_timestamp, mga_sunrise_timestamp, mga_sunset_timestamp)
        minutes_in_mga_temporal_hour = str(round(seconds_in_mga_temporal_hour / 60,1)) # str(convert_seconds(seconds_in_mga_temporal_hour))
    else:
        
        mga_temporal_time = reverse("שגיאה  ")
        minutes_in_mga_temporal_hour = 0.0


    ###################### חישוב שלב הירח הנוכחי #####################
    
    MoonPhase.tim = round(current_timestamp) ############### אם לא מגדירים את זה אז הזמן הוא לפי הזמן הפנימי של הבקר
    mp = MoonPhase()  # במיקרופייתון צריך לעשות mp = MoonPhase(lto=location_offset_hours) אך בפייתון רגיל צריך לעשות בלי lto כדי שיהיה מדוייק בדומה לסקייפילד ואינני יודע כעת מדוע
    phase = mp.phase()
    phase_percent = round(phase * 100,1)
                      
    ##################################################################################################################3
    # משתנה שמחזיר טרו אם הזמן הוא כרגע אחרי השקיעה ולפני 12 בלילה ולכן התאריך העברי הנכון הוא התאריך הלועזי הבא
    heb_date_is_next_greg_date = sunset and current_timestamp > sunset and current_timestamp > sunrise # current_timestamp > sunrise אומר שמדובר לפני 12 בלילה
    
    # אם התאריך העברי מקביל לתאריך הלועזי של מחר כי מדובר אחרי השקיעה ולפני 12 בלילה מחשבים את הנתונים על מחר
    if heb_date_is_next_greg_date:    
        # חישוב התאריך הלועזי של מחר כלומר בדיוק עוד 24 שעות. זה נדרש כי התאריך העברי מהשקיעה עד 12 בלילה שווה לתאריך הלועזי של מחר
        tomorrow_tm = time.gmtime(current_timestamp+86400) # יש 86400 שניות ביממה
        g_year, g_month, g_day, g_rtc_week_day = (tomorrow_tm[0], tomorrow_tm[1], tomorrow_tm[2], tomorrow_tm[6])
    
    # בכל מקרה אחר התאריך העברי מקביל לתאריך הלועזי הנוכחי
    else:
        g_year, g_month, g_day, g_rtc_week_day = year, month, day, rtc_week_day
    
    
    # חישוב תאריך עברי נוכחי באמצעות ספרייה ייעודית. כמו כן מחושב האם מדובר בחג
    heb_date_string, tuple_heb_date, holiday_name, lite_holiday_name, is_rosh_chodesh = get_heb_date_and_holiday_from_greg_date(g_year, g_month, g_day)
    # חישוב היום בשבוע המתאים לתאריך העברי הנכון לרגע זה
    heb_weekday = get_normal_weekday(g_rtc_week_day)
    # שם בעברית של היום העברי המתאים לתאריך העברי הנוכחי המוגדר משקיעה מישורית לשקיעה מישורית
    heb_weekday_string = heb_weekday_names(heb_weekday)

    ##############################################################################
    # איזור שאחראי להגדיר ששעון ההלכה לא ייכנס אוטומטית למצב שינה בשבת ובחג. אך לא מחושב יום טוב שני
    
    # חישוב האם שבת. שבת מוגדרת מהשקיעה של סוף יום שישי עד השקיעה של סוף שבת
    is_shabat = heb_weekday == 7
    
    # חישוב תוספות לשבת כלומר מיום שישי חצי שעה לפני השקיעה עד השקיעה וכן בשבת מהשקיעה ועד צאת שבת שבלוחות
    normal_weekday = get_normal_weekday(rtc_week_day) # חישוב היום בשבוע של התאריך הלועזי בדווקא
    half_hour_before_sunset_until_sunset =  sunset and current_timestamp >= (sunset - 1800) and current_timestamp < sunset # 1800 שניות זה חצי שעה לפני השקיעה
    sunset_until_motsaei_shabat_luchot = sunset and current_timestamp > sunset and s_alt > -8.5
    is_tosafot_leshabat = (normal_weekday == 6 and half_hour_before_sunset_until_sunset) or (normal_weekday == 7 and sunset_until_motsaei_shabat_luchot)
    
    ##############################################################################
      
    # מכאן והלאה ההדפסות למסך
    
    #voltage_string = f"{round(voltage,1)}v"
    voltage_string = f"**%"
    
    greg_date_string = f'{day:02d}/{month:02d}/{year:04d}' 
    time_string = f'{hour:02d}:{minute:02d}:{second:02d}'
    
    # מהשקיעה עד 12 בלילה מוסיפים את המילה ליל כי היום בשבוע והתאריך העברי מקבלים לתאריך הלועזי של מחר
    leil_string = reverse("ליל: ") if heb_date_is_next_greg_date else ""
    # אם אין שעון והוגדר זמן שרירותי או שהשעה נלקחה מהשעון הפנימי שכנראה אינו מדוייק מוסיפים סימני קריאה אחרי התאריך העברי
    heb_date_to_print = f'{leil_string}{reverse(heb_weekday_string)}, {reverse(heb_date_string)}'
    # בלינוקס צריך לשנות את הסדר בגלל בעיית תצוגת עברית
    if not is_windows:
        heb_date_to_print = f'{reverse(heb_date_string)} ,{reverse(heb_weekday_string)}{leil_string}'
    #magrab_time = calculate_magrab_time(current_timestamp, sunset_timestamp) if sunrise else reverse("שגיאה  ") # רק אם יש זריחה ושקיעה אפשר לחשב
    utc_offset_string = 'utc+00' if location_offset_hours == 0 else f'utc+{location_offset_hours:02}' if location_offset_hours >0 else f'utc-{abs(location_offset_hours):02}'
    coteret = f'{voltage_string} - {reverse("שעון ההלכה")} - {reverse(location["heb_name"])}*'
    # בלינוקס צריך לשנות את הסדר בגלל בעיית תצוגת עברית
    if not is_windows:
        coteret = f'{voltage_string} - {reverse(location["heb_name"])} - {reverse("שעון ההלכה")}*'
    
    
    # עדכון שורת הכותרת
    canvas.itemconfig(title_id, text=coteret)
    
    # עדכון התאריך העברי כולל שינוי צבע הרקע של התאריך העברי לפי המשתנה
    
    # איזור תאריך עברי כולל צבע מתאים לימי חול ולשבתות וחגים
    # צבע הטקסט והרקע של התאריך העברי: ביום חול לבן על שחור ובשבת וחג שחור על צהוב, ובחגים דרבנן כולל תעניות שחור על ציאן
    
    #################################################################################
    # איזור לטיפול בכך שהרקע הצבעוני שתחת התאריך העברי יהיה רק תחתיו ולא בכל השורה
    
    # קבלת טקסט עם הנתונים על הפונט שהוגדר לאיזור התאריך העברי והפרדתו למערך של שלושה נתונים
    heb_date_font_tuple = canvas.itemcget(heb_date_id, "font").split()
    # הפרדת שלושת הנתונים
    font_family, font_size, font_weight = heb_date_font_tuple
    # בניית אובייקט פונט של טקינטר מתוך מערך שלושת הנתונים
    heb_date_font = font.Font(family=font_family, size=int(font_size), weight=font_weight)   
    # מדידת אורך הטקסט של התאריך העברי בהתאם לפונט שבו השתמשנו לכתוב אותו
    text_width = heb_date_font.measure(heb_date_to_print)
    y1 = 21 * scale # גובה התחלת הרקע מראש המסך
    y2 = 38 * scale # גובה סיום הרקע מראש המסך
    center_x = 160 * scale # מרכז המסך מימין לשמאל
    padding = 30 # שוליים נוספים מעבר לטסקט אם רוצים שהרקע יחרוג קצת. לא הכרחי אפשר לעשות 0
    x1 = center_x - (text_width / 2) - padding # מרחק תחילת הטקסט ממרכז המסך
    x2 = center_x + (text_width / 2) + padding # מרחק סוף הטקסט ממרכז המסך
    canvas.coords(heb_date_rect_id, x1, y1, x2, y2) # עדכון הרקע כך שיהיה רק בגבולות של הטקסט

    ####################################################################
    

    HEB_DATE_FG, HEB_DATE_BG  = ("black", "yellow") if is_shabat or holiday_name else ("black", "cyan") if lite_holiday_name or is_rosh_chodesh else ("white", "black")
    canvas.itemconfig(heb_date_rect_id, fill=HEB_DATE_BG)
    canvas.itemconfig(heb_date_id, text=heb_date_to_print, fill=HEB_DATE_FG)
    
    # עדכון שעה זמנית גרא ומגא ודקות בשעה זמנית לכל אחת מהשיטות
    canvas.itemconfig(minutes_in_gra_temporal_hour_id, text=minutes_in_temporal_hour)
    canvas.itemconfig(gra_temporal_hour_id, text=temporal_time)
    canvas.itemconfig(minutes_in_mga_temporal_hour_id, text=minutes_in_mga_temporal_hour)
    canvas.itemconfig(mga_temporal_hour_id, text=mga_temporal_time)
    
    # עדכון מידע על השמש והירח גובה אזימוט ושלב בחודש
    canvas.itemconfig(sun_az_id, text=f'{"  " if s_az < 10 else " " if s_az < 100 else ""}{round(s_az)}°')
    canvas.itemconfig(sun_alt_id, text=f'{" " if s_alt > 0 else ""}{" " if abs(s_alt) <10 else ""}{round(s_alt,3):.3f}°') 
    canvas.itemconfig(moon_az_id, text=f'{"  " if m_az < 10 else " " if m_az < 100 else ""}{round(m_az)}°')
    canvas.itemconfig(moon_alt_id, text=f' {" " if m_alt > 0 else ""}{" " if abs(m_alt) <10 else ""}{m_alt:.3f}°')
    canvas.itemconfig(moon_phase_id, text=f'    {phase_percent:.1f}%')
    
    # עדכון שורת ההסברים
    global current_screen_halach_clock 
    text = reverse(esberim[int(current_screen_halach_clock)][0])  # רוורס של הטקסט העברי
    time_value = esberim[int(current_screen_halach_clock)][1]  # הערך להצגה
    CCC = f"{time_value}  :{text}" if time_value != "" else f"{text}"
    canvas.itemconfig(hesberim_id, text=f"{CCC}")
    current_screen_halach_clock = (current_screen_halach_clock + 0.25) % len(esberim)  # זה גורם מחזור של שניות לאיזה נתונים יוצגו במסך
    
    
    ############################################################################
    ###########################################################################3
    ############################################################################
    # איזור הדפסת זמנים בשעון רגיל. כל האיזור עדיין בבניה
    
    #חישוב מספר השקיעה מהזריחה לשקיעה
    seconds_day_gra = (sunset - sunrise) / 12 if sunrise and sunset else None
    seconds_day_mga = (mga_sunset - mga_sunrise) / 12 if mga_sunrise and mga_sunset else None
    
    def hhh(start_time, seconsd_per_hour, hour):
        if seconsd_per_hour:
            AAA = start_time + (seconsd_per_hour * hour)
             # עיגול לדקה הקרובה
            total_seconds = int(AAA + 30) // 60 * 60 
            time_value = time.gmtime(total_seconds)
            return time.strftime("%H:%M", time_value)
            # אם רוצים בלי עיגול אלא כולל שניות
            #time_value = time.gmtime(AAA)
            #return time.strftime("%H:%M:%S", time_value) 
        else:
            return reverse("שגיאה")
        
    

    zmanim = [
        
        [f"עלות השחר {reverse(16)}", hhh(mga_sunrise, seconds_day_mga, hour=0)],
        ["זריחה", hhh(sunrise, seconds_day_mga, hour=0)],
        ["סוף שמע מגא", hhh(mga_sunrise, seconds_day_mga, hour=3)],
        ["סוף שמע גרא", hhh(sunrise, seconds_day_gra, hour=3)],
        ["סוף תפילה מגא",  hhh(mga_sunrise, seconds_day_mga, hour=4)],
        ["סוף תפילה גרא", hhh(sunrise, seconds_day_gra, hour=4)],
        ["חצות", hhh(sunrise, seconds_day_gra, hour=6)],
        ["מנחה גדולה", hhh(sunrise, seconds_day_gra, hour=6.5)],
        ["מנחה קטנה", hhh(sunrise, seconds_day_gra, hour=9.5)],
        ["פלג המנחה", hhh(sunrise, seconds_day_gra, hour=10.75)],
        ["שקיעה", hhh(sunrise, seconds_day_gra, hour=12)],
        [f"צאת הכוכבים דרבינו תם {reverse(16)}", hhh(mga_sunrise, seconds_day_mga, hour=12)],
    ]

    global current_screen_zmanim
    # אם רוצים זמן אחד בשורה
    #text = reverse(zmanim[int(current_screen_zmanim)][0])
    #time_i = zmanim[int(current_screen_zmanim)][1]
    #SSS = f' {text}: {time_i}'
    #canvas.itemconfig(sgb_id, text=SSS)
    #current_screen_zmanim = (current_screen_zmanim + 0.25) % len(zmanim)
    
    # אם רוצים שני זמנים בשורה
    lines = []
    base_index = int(current_screen_zmanim) * 2

    for i in range(2):
        index = (base_index + i) % len(zmanim)
        label = reverse(zmanim[index][0])
        time_val = zmanim[index][1]
        # בלינוקס צריך סדר הפוך מווינדוס בגלל בעיית הצגת עברית
        if is_windows:
            lines.append(f'{label}: {time_val}')
        else:
            lines.append(f'{time_val} :{label}')

    SSS = '   |   '.join(lines)
    canvas.itemconfig(zmanim_id, text=SSS)
    current_screen_zmanim = (current_screen_zmanim + 0.15) % ((len(zmanim) + 1) // 2)  # חלוקה לשלשות, מעוגלת כלפי מעלה

    
    #############################################################################
    #############################################################################
    #############################################################################
    
    # עדכון תאריך לועזי שעה רגילה ואיזור זמן
    canvas.itemconfig(utc_offset_id, text=utc_offset_string)
    canvas.itemconfig(time_id, text=time_string)
    canvas.itemconfig(greg_date_id, text=greg_date_string)
    
    gc.collect() # ניקוי הזיכרון חשוב נורא כדי למנוע קריסות
    
    ####################################################################
    # כשהמונה מגיע ל- 45 מדמים לחיצה על מקש שיפט במקלדת כדי למנוע כיבוי אוטומטי  
    global counter_shift
    if counter_shift >= 45:
        #print("counter_shift_press")
        keyboard = keyboard_Controller()
        keyboard.press(keyboard_Key.shift)
        keyboard.release(keyboard_Key.shift)
        counter_shift = 0.0 # איפוס המונה
    
    counter_shift += 1 # בכל שנייה המונה מתקדם באחד
    
    ####################################################################
               
    # חזרה על העדכון כל שנייה מחדש
    root_hw.after(1000, main_halach_clock)
    
    
######################################################################################################33

# לולאת רענון חשובה ביותר שחוזרת על עצמה כל הזמן והיא זו שמפעילה את הפונקצייה הראשית כל שנייה מחדש
# התחלת העדכון
main_halach_clock()

# בודק האם מותקן במחשב גופן מירים שדרוש לתצוגה טובה של שעון ההלכה
is_miriam = "Miriam" in font.families()
if not is_miriam:
    messagebox.showinfo(reverse("התרעה"), reverse("הגופן העברי 'מירים' מווינדוס חייב להיות מותקן לצורך תצוגה נכונה"))
    root_hw.focus_force()  # מחזיר את הפוקוס לחלון הראשי כי אחרת אי אפשר לצאת מהתוכנה בלחיצה על אסקייפ או לשנות מיקומים בלחצני החיצים

# הפעלת החלון הראשי בקביעות
root_hw.mainloop()



