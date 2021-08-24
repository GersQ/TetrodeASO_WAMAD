LIBRARY IEEE;
USE IEEE.STD_LOGIC_1164.ALL;
USE IEEE.NUMERIC_STD.ALL;
USE IEEE.MATH_REAL.LOG2;
USE IEEE.MATH_REAL.CEIL;



USE WORK.VHDLTOOL.ALL;

ENTITY SMOOTH_CONTROLLER IS 
GENERIC (
				K:INTEGER; 
				CHANNEL: INTEGER
			);
PORT(
		CLK,RESET,ENABLE	: IN STD_LOGIC;
		READY			: OUT STD_LOGIC;
		INDEXW		 	: OUT UNSIGNED(INTEGER(CEIL(LOG2(REAL(4*K*CHANNEL+1))))-1 DOWNTO 0);
		INDEXR			: OUT POINTERS (0 TO 4*K)
		);
END ENTITY SMOOTH_CONTROLLER;

ARCHITECTURE RTL OF SMOOTH_CONTROLLER IS
SIGNAL INDEXW_TMP	:  UNSIGNED(INTEGER(CEIL(LOG2(REAL(4*K*CHANNEL+1))))-1 DOWNTO 0);
SIGNAL FLAG		: STD_LOGIC;

BEGIN

WRITING_MEM:PROCESS(CLK,RESET)
VARIABLE CNT: UNSIGNED(INTEGER(CEIL(LOG2(REAL(4*K*CHANNEL+1))))-1 DOWNTO 0);
BEGIN
IF(RESET='0') THEN
	CNT:=(OTHERS=>'0');
	FLAG<='1';
ELSIF(RISING_EDGE(CLK)) THEN

IF(ENABLE='0') THEN
	CNT:=CNT+1;
	IF(CNT=(4*K*CHANNEL+1)) THEN
		FLAG<='0';
			CNT:=(OTHERS=>'0');
			END IF;
END IF;
END IF;
INDEXW_TMP<=CNT;
END PROCESS WRITING_MEM;

READY<=FLAG OR ENABLE;
INDEXW<=INDEXW_TMP;

-- This process creates the 4*k+1 read address.
READING_MEM:PROCESS(CLK,RESET)
VARIABLE TMP: UNSIGNED(INTEGER(CEIL(LOG2(REAL(4*K*CHANNEL+1))))-1 DOWNTO 0);
VARIABLE R: POINTERS (0 TO 4*K);
BEGIN
IF(RESET='0') THEN

		
	FOR I IN 0 TO 4*K LOOP
	
		TMP:=TO_UNSIGNED(I*CHANNEL,INDEXW'LENGTH);
		R(I):=TMP;
		
	END LOOP;
ELSIF(RISING_EDGE(CLK)) THEN

IF(FLAG='0') THEN
	

 	
			
	FOR I IN 0 TO 4*K LOOP
		R(I):=R(I)+1;

		IF(R(I)=(4*K*CHANNEL+1)) THEN
			R(I):=(OTHERS=>'0');
		END IF;
	END LOOP;
			
	
	
			
END IF;
END IF;
INDEXR<=R;
END PROCESS READING_MEM;
END ARCHITECTURE RTL ;
