LIBRARY IEEE;
USE IEEE.STD_LOGIC_1164.ALL;
USE IEEE.NUMERIC_STD.ALL;
USE IEEE.MATH_REAL.CEIL;
USE IEEE.MATH_REAL.LOG2;

--SIMULANEOUS R/W OPERATION TO REDUCE AS MUCH AS POSSIBLE THE LOGIC REQUIRED

ENTITY BUTT_CONTROLLER IS
GENERIC (CHANNEL: INTEGER:=1);
PORT(	
		CLK,RESET, ENABLE	:  IN STD_LOGIC;
		READY					:	OUT STD_LOGIC;
		INDEXW				:  OUT UNSIGNED(INTEGER(CEIL(LOG2(REAL(CHANNEL*2))))-1 DOWNTO 0);
		INDEXR_OLD			:  OUT UNSIGNED(INTEGER(CEIL(LOG2(REAL(CHANNEL*2))))-1 DOWNTO 0)
		);
		
END ENTITY BUTT_CONTROLLER;


ARCHITECTURE RTL OF BUTT_CONTROLLER IS
SIGNAL FULL,READY_TMP: STD_LOGIC;



BEGIN
 
INDEXING_MEMWRITE:PROCESS(CLK,RESET)
VARIABLE CNT: UNSIGNED(INTEGER(CEIL(LOG2(REAL(CHANNEL*2))))-1 DOWNTO 0);
BEGIN
IF (RESET='0') THEN
	CNT:=(OTHERS=>'0');
	FULL<='1';
ELSIF (RISING_EDGE(CLK)) THEN
	IF(ENABLE='0') THEN
		CNT:=CNT+1;
		
		IF(CNT=CHANNEL*2-1)THEN
			FULL<='0';
		END IF;
		
	END IF;
END IF;
INDEXW<=CNT;
END PROCESS INDEXING_MEMWRITE;

PROCESS(CLK,RESET)
BEGIN
IF(RESET='0') THEN
READY_TMP<='1';
ELSIF(RISING_EDGE(CLK)) THEN
READY_TMP<=FULL;
END IF;
END PROCESS;

--IF A COMBINATORIAL OUTPUT IS AN ISSUE THEN DELAY THE READY_TMP AND SET IT LOW WHEN THE OUTPUT OF REG IS EVALUATED.
READY<=READY_TMP;


INDEXING_READMEMOLD:PROCESS(CLK,RESET)
VARIABLE CNT: UNSIGNED(INTEGER(CEIL(LOG2(REAL(CHANNEL*2))))-1 DOWNTO 0);
BEGIN
IF (RESET='0') THEN
	CNT:=TO_UNSIGNED(CHANNEL-1,CNT'LENGTH);

ELSIF (RISING_EDGE(CLK)) THEN
	IF(FULL='0') THEN
		CNT:=CNT+1;
		
		
END IF;
END IF;
INDEXR_OLD<=CNT;
END PROCESS INDEXING_READMEMOLD;

END ARCHITECTURE RTL;
