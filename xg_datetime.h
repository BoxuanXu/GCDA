/*
 * $Id$
 */
#ifndef XG_DATETIME_H
#  define XG_DATETIME_H

#include <time.h>

int monthlength[]={
    -1,
    31,28,31,30,31,30,31,31,30,31,30,31
};
const char* monStr[] = {"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"};

bool isleapyear(int y){
    if(y % 4 == 0){
        if(y % 100 == 0)
            return (y % 400 ==0);
        else 
            return true;
    }
    return false;
}

clock_t Cal_time(clock_t& start_time)
{
	clock_t time_now = clock();
	cout<< "Running time is: "<<static_cast<double>(time_now-start_time)/CLOCKS_PER_SEC*1000<<"ms"<<endl;//print during time
	
	time_t timep;
	time(&timep);
	printf("Now time is : %s\n", ctime(&timep));

	time_now = clock();	
	return time_now;
}	

struct datetime{
    int year, month, day, hour, minute, seconds;
    /*
     * year     -> YYYY, YY
     * month    -> MM, M
     * day      -> DD, D
     * hour     -> hh, h
     * minute   -> mm, m
     * seconds  -> ss, s
     */
    double msecond;
    datetime(int y = 1970, int M = 1, int d = 1, int h=0, int m=0, int s=0):
        year(y), month(M), day(d), hour(h), minute(m), seconds(s){
    }
    datetime(const string& datestr, const string& frmt = "YYYY-MM-DD hh:mm:ss"){
        year = 0; month = 0; day = 0; hour = 0; minute = 0; seconds = 0;
        assert(datestr.size() == frmt.size());
        int s = datestr.size(), i = 0;
        while(i < s){
            if(i + 4 <= s && frmt.substr(i,4) == "YYYY"){
                year = atoi(datestr.substr(i, 4).c_str());
                i += 4;
            }
            else if(i + 2 <= s && frmt.substr(i,2) == "YY"){
                year = atoi(datestr.substr(i, 2).c_str());
                year += year > 50 ? 1900 : 2000;
                i += 2;
            }
            else if(i + 2 <= s && frmt.substr(i,2) == "MM"){
                month = atoi(datestr.substr(i, 2).c_str());
                i += 2;
            }
            else if(i + 2 <= s && frmt.substr(i,2) == "DD"){
                day = atoi(datestr.substr(i, 2).c_str());
                i += 2;
            }
            else if(i + 2 <= s && frmt.substr(i,2) == "hh"){
                hour = atoi(datestr.substr(i, 2).c_str());
                i += 2;
            }
            else if(i + 2 <= s && frmt.substr(i,2) == "mm"){
                minute = atoi(datestr.substr(i, 2).c_str());
                i += 2;
            }
            else if(i + 2 <= s && frmt.substr(i,2) == "ss"){
                seconds = atoi(datestr.substr(i, 2).c_str());
                //cout << "seconds = " << seconds  << endl;
                i += 2;
            }
            else if(i + 1 <= s && frmt.substr(i,1) == "M"){
                month = atoi(datestr.substr(i, 1).c_str());
                i += 1;
            }
            else if(i + 1 <= s && frmt.substr(i,1) == "D"){
                day = atoi(datestr.substr(i, 1).c_str());
                i += 1;
            }
            else if(i + 1 <= s && frmt.substr(i,1) == "h"){
                hour = atoi(datestr.substr(i, 1).c_str());
                i += 1;
            }
            else if(i + 1 <= s && frmt.substr(i,1) == "m"){
                minute = atoi(datestr.substr(i, 1).c_str());
                i += 1;
            }
            else if(i + 1 <= s && frmt.substr(i,1) == "s"){
                seconds = atoi(datestr.substr(i, 1).c_str());
                //cout << "seconds = " << seconds  << endl;
                i += 1;
            }
            else{
                i += 1;
            }
        }
    }
    /*
     * unit: day
     * since:12h Jan 1, 4713 BC
     */
    double JD()const{//unit : day
        int a = (14 - month) / 12;
        int y = year + 4800 - a;
        int m = month + 12 * a - 3;
        return day + (153 * m + 2 ) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045 + (hour - 12) / 24.0 + minute / 1440.0 + seconds / 86400.0;
    }
    /*
     * unit: seconds
     * since:0h Jan 1, 1970
     */
    long long unixtime()const{ 
        int a = (14 - month) / 12;
        int y = year + 4800 - a;
        int m = month + 12 * a - 3;
        return (day + (153 * m + 2 ) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045 - 2440588) * 86400 + hour * 3600 + minute * 60 + seconds ;
    }
    /*
     * TAU value (# of hours since 1 Jan 1985)
     */
    double tau()const{
        return 24 * (JD() - datetime(1985,1,1).JD());
    }
    int monthmax()const{
        return monthlength[month] + ((month == 2) ? isleapyear(year):0);
    }
    void uniform(){
        /* second */
        while(seconds >= 60){
            seconds -= 60;
            minute += 1;
        }
        while(seconds < 0){
            seconds += 60;
            minute -= 1;
        }
        /* minute */
        while(minute >= 60){
            minute -= 60;
            hour += 1;
        }
        while(minute < 0){
            minute += 60;
            hour -= 1;
        }
        /* hour */
        while(hour >= 24){
            hour -= 24;
            day += 1;
        }
        while(hour < 0){
            hour += 24;
            day -= 1;
        }
        /*****************************/
        /* day                       */
        while(true){
            int l = monthlength[month] + ((month == 2) ? isleapyear(year):0);
            //printf("l = %d\n", l );
            if(day > l){
                day -= l;
                month += 1;
                while(month > 12){
                    month -= 12;
                    year += 1;
                }
            }
            else 
                break;
        }
        /* month */
        
    }
    string str(const string& frmt = "YYYY-MM-DD hh:mm:ss")const{
        string ans = "";
        int i = 0, s = frmt.size();
        char sbuf[5];
        while(i < s){
            if(i + 4 <= s && frmt.substr(i,4) == "YYYY"){
                sprintf(sbuf, "%.4d", year);
                i += 4;
            }
            else if(i + 2 <= s && frmt.substr(i,2) == "YY"){
                sprintf(sbuf, "%.2d", year - year / 100 * 100);
                i += 2;
            }
            else if(i + 2 <= s && frmt.substr(i,2) == "MM"){
                sprintf(sbuf, "%.2d", month);
                i += 2;
            }
            else if(i + 2 <= s && frmt.substr(i,2) == "DD"){
                sprintf(sbuf, "%.2d", day);
                i += 2;
            }
            else if(i + 2 <= s && frmt.substr(i,2) == "hh"){
                sprintf(sbuf, "%.2d", hour);
                i += 2;
            }
            else if(i + 2 <= s && frmt.substr(i,2) == "mm"){
                sprintf(sbuf, "%.2d", minute);
                i += 2;
            }
            else if(i + 2 <= s && frmt.substr(i,2) == "ss"){
                sprintf(sbuf, "%.2d", seconds);
                i += 2;
            }
            else if(i + 1 <= s && frmt.substr(i,1) == "M"){
                sprintf(sbuf, "%d", month);
                i += 1;
            }
            else if(i + 1 <= s && frmt.substr(i,1) == "D"){
                sprintf(sbuf, "%d", day);
                i += 1;
            }
            else if(i + 1 <= s && frmt.substr(i,1) == "h"){
                sprintf(sbuf, "%d", hour);
                i += 1;
            }
            else if(i + 1 <= s && frmt.substr(i,1) == "m"){
                sprintf(sbuf, "%d", minute);
                i += 1;
            }
            else if(i + 1 <= s && frmt.substr(i,1) == "s"){
                sprintf(sbuf, "%d", seconds);
                i += 1;
            }
            else{
                sprintf(sbuf, "%c", frmt[i]);
                i += 1;
            }
            ans += sbuf;
        }
        return ans;
    }
    string YYYYMMDDhhmmss(){
        char buf[15];
        sprintf(buf, "%.4d%.2d%.2d%.2d%.2d%.2d", year, month, day, hour, minute, seconds);
        return buf;
    }
    string YYYYMMDDhhmm(){
        char buf[11];
        sprintf(buf, "%.4d%.2d%.2d%.2d%.2d", year, month, day, hour,minute);
        return buf;
    }
    string YYYYMMDD(){
        char buf[9];
        sprintf(buf, "%.4d%.2d%.2d", year, month, day);
        return buf;
    }
};


struct timespan{
    int year, month, day, hour, minute, seconds;
    double msecond;
    timespan(int y = 0, int M = 0, int d = 0, int h=0, int m=0, int s=0):
        year(y), month(M), day(d), hour(h), minute(m), seconds(s){
    }
};

/*
timespan operator-(const datetime& d1, const datetime& d2){
}
*/

datetime operator-(datetime d1, const timespan& s){
    d1.year    -= s.year;
    d1.month   -= s.month;
    d1.day     -= s.day;
    d1.hour    -= s.hour;
    d1.minute  -= s.minute;
    d1.seconds -= s.seconds;
    d1.uniform();
    return d1;
}

datetime operator+(datetime d1, const timespan& s){
    d1.year    += s.year;
    d1.month   += s.month;
    d1.day     += s.day;
    d1.hour    += s.hour;
    d1.minute  += s.minute;
    d1.seconds += s.seconds;
    d1.uniform();
    return d1;
}
datetime operator+=(datetime& d1, const timespan& s){
    d1 = d1 + s;
    return d1;
}
bool operator==(const datetime & d1, const datetime & d2){
	if(d1.year == d2.year && d1.month == d2.month && d1.day == d2.month && d1.hour == d2.hour && d1.minute == d2.minute && d1.seconds == d2.seconds)
		return true;
	else
		return false;
}

bool operator<(const datetime & d1, const datetime & d2){
    if(d1.year != d2.year) return d1.year < d2.year;
    if(d1.month != d2.month) return d1.month < d2.month;
    if(d1.day != d2.day) return d1.day < d2.day;
    if(d1.hour != d2.hour) return d1.hour < d2.hour;
    if(d1.minute != d2.minute) return d1.minute < d2.minute;
    return d1.seconds < d2.seconds;
}
bool operator>(const datetime & d1, const datetime & d2){
    return d2 < d1;
}
bool operator<=(const datetime & d1, const datetime & d2){
    return !(d1 > d2);
}
bool operator>=(const datetime & d1, const datetime & d2){
    return !(d1 < d2);
}


#endif /* ifndef XG_DATETIME_H */

