scalar
   t 'number of roaming conductors' /4/
   we 'number of working events to visit' /5/
   v 'roaming conductor speed' /40/
   h 'shift hours' /12/
;
sets
   i working events /i1*i5/
   c roaming conductors /c1*c3/
;
alias (i,j);


display j;

Parameters
       a(i)  time window start
         /    i1   0
              i2   0
              i3   0
              i4   0
              i5   0  /

       b(i)  time window end
         /    i1   12
              i2   12
              i3   12
              i4   12
              i5   12  /
;


Table d(i,j) 'distance matrix (KM)'
      i1   i2   i3   i4   i5 
 i1    0  137  125   39    0 
 i2  137    0  129  133  137   
 i3  125  129    0  229  125  
 i4   39  133  229    0   39   
 i5    0  137  125   39    0  

;



Positive variables w(i,c);

Binary variables  x(i,j,c);

variables z;
Equations
depotstart,depotend,inout,shifthour,visitonce1,visitonce2,timewindow1,timewindow2,sequence,sub1,sub2,sub3,sub4,sub5,workhour1,workhour2,depart,return,num;


depotstart(c).. sum(j, x('i1',j,c)) =l= 1;
depotend(c)..   sum(j, x(j,'i5',c)) =l= 1;

inout(i,c)$(ord(i)<>1 and ord(i)<>5).. sum(j$(ord(j)<>5),x(j,i,c))- sum(j$(ord(j)<>1),x(i,j,c)) =e= 0;

shifthour(c).. sum((i,j),x(i,j,c)*d(i,j)/v) =l= h;

visitonce1(i)$(ord(i)<>1 and ord(i)<>5).. sum((j,c)$(ord(j)<>1),x(i,j,c)) =e= 1;
visitonce2(j)$(ord(j)<>1 and ord(j)<>5).. sum((i,c)$(ord(i)<>5),x(i,j,c)) =e= 1;

timewindow1(i,c)$(ord(i)<>1 and ord(i)<>5).. w(i,c) =g= a(i)*sum(j,x(i,j,c));
timewindow2(i,c)$(ord(i)<>1 and ord(i)<>5).. w(i,c) =l= b(i)*sum(j,x(i,j,c));

sequence(i,j,c)$(ord(i)<>ord(j)).. w(i,c) + d(i,j)/v - w(j,c) =l= 17.725*(1-x(i,j,c));


sub1(c)..x('i1','i1',c)=e=0;
sub2(c)..x('i2','i2',c)=e=0;
sub3(c)..x('i3','i3',c)=e=0;
sub4(c)..x('i4','i4',c)=e=0;
sub5(c)..x('i5','i5',c)=e=0;

workhour1(c)..w('i1',c) =g= 0;
workhour2(c)..w('i5',c) =l= h;

depart(j,c).. x('i5',j,c) =e= 0;
return(i,c).. x(i,'i1',c) =e= 0;

num.. z  =e=  sum((c,j), x('i1',j,c));


Model conductor /all/ ;

Solve conductor using miP minimizing z ;
