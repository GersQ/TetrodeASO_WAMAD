LIBRARY IEEE;
USE IEEE.STD_LOGIC_1164.ALL;
USE IEEE.NUMERIC_STD.ALL;
USE IEEE.MATH_REAL.ALL;


--3 POINTERS TO WRITE AND READ, SIMULATANEOUS W/R OPERATION

ENTITY BUTTERWORTH_BUFFV IS

GENERIC(
			CHANNEL	: INTEGER :=30; 
			N			: INTEGER :=10
			); 
			
PORT(
		CLK,RESET			: IN STD_LOGIC;
		READY					: IN STD_LOGIC;
		SAMPLE				: IN STD_LOGIC_VECTOR(N-1 DOWNTO 0);    --QN-1.0
		FILTERED				: OUT STD_LOGIC_VECTOR(N-1 DOWNTO 0); --QN-1.0

		INDEXW,INDEXR_OLD	: IN	UNSIGNED(INTEGER(CEIL(LOG2(REAL(CHANNEL*2))))-1 DOWNTO 0)

		); 
END ENTITY BUTTERWORTH_BUFFV;


ARCHITECTURE RTL OF BUTTERWORTH_BUFFV IS
TYPE MAT1 IS ARRAY(0 TO 2*CHANNEL-1) OF STD_LOGIC_VECTOR(N-1 DOWNTO 0);
SIGNAL Q,S: MAT1;
SIGNAL SUM,M1,M3,M4,M5,ADD1,ADD2: SIGNED (2*N-1 DOWNTO 0);
SIGNAL X: STD_LOGIC_VECTOR(N-1 DOWNTO 0);



--FILTER CONSTANTS IIR BUTTERWORTH II ORDER PASSBAND 300HZ-3000HZ.
--H(Z)=(BO+B1*Z^-1+B2*Z^-2) / (AO+A1*Z^-1+A2*Z^-2)
-- Filter Coefficient  sampling rate, 20kHz Q1.8
CONSTANT B0: SIGNED(N-1 DOWNTO 0):= TO_SIGNED(49,10); --Q9.0
CONSTANT B2: SIGNED(N-1 DOWNTO 0):= TO_SIGNED(-49,10); --Q9.0
CONSTANT A1: SIGNED(N-1 DOWNTO 0):= TO_SIGNED(-336,10); --Q9.0
CONSTANT A2: SIGNED(N-1 DOWNTO 0):= TO_SIGNED(97,10); --Q9.0

BEGIN

--PROCESS TO SAMPLE THE INPUT 
START:PROCESS(CLK,RESET)
BEGIN
IF (RESET='0') THEN
	X<=(OTHERS=>'0');

ELSIF (RISING_EDGE(CLK)) THEN
	X<=SAMPLE;

END IF;

END PROCESS START;


BUFFIN:PROCESS(CLK,RESET)
BEGIN
IF (RESET='0') THEN
	Q<=(OTHERS=>(OTHERS=>'0'));
ELSIF (RISING_EDGE(CLK)) THEN
	Q(TO_INTEGER(INDEXW))<=X;
END IF;
END PROCESS BUFFIN;


-- MULTIPLICATIONS
M1<=SIGNED(X)*B0;
M3<=SIGNED(Q(TO_INTEGER(INDEXW)))*B2;
M4<=SIGNED(S(TO_INTEGER(INDEXR_OLD)))*A1;
M5<=SIGNED(S(TO_INTEGER(INDEXW)))*A2;

ADD1<=M1+M3;
ADD2<=M4+M5;
SUM<=ADD1-ADD2; -- OR CHANGE THE SIGN OF A COEFFICIENT 
FILTERED<=STD_LOGIC_VECTOR(SUM(2*N-3 DOWNTO N-2)) WHEN READY='0' ELSE (OTHERS=>'0'); --Q19.0->Q10.0 SHIFT RX 9 POSITION THEN EXTRACT 10 RIGHTMOST LSB TO GET Q9.0

BUFFOUT:PROCESS(CLK,RESET)
BEGIN
IF (RESET='0') THEN
	S<=(OTHERS=>(OTHERS=>'0'));
ELSIF (RISING_EDGE(CLK)) THEN
	IF(READY='0') THEN  --POTREI PURE ELIMINARLO
	 S(TO_INTEGER(INDEXW))<=STD_LOGIC_VECTOR(SUM(2*N-3 DOWNTO N-2));
	END IF;
END IF;
END PROCESS BUFFOUT;
	


END ARCHITECTURE RTL;
