function TotalSeconds = timeToSec(hours, mins, secs)
hoursSec = hours * 60*60;
minsSec = mins *60;
TotalSeconds = hoursSec + minsSec + secs;