LIBRARY IEEE;
USE IEEE.STD_LOGIC_1164.ALL;
USE IEEE.NUMERIC_STD.ALL;
USE IEEE.MATH_REAL.LOG2;
USE IEEE.MATH_REAL.CEIL;

--LATEST VERSION INCLUDING A BETTER APPROACH FOR AA, WITH


ENTITY WA IS 
	GENERIC (
					WINDOW: INTEGER; 
					NBIT: INTEGER
                                        				);
	PORT(	
			SAMPLE: IN STD_LOGIC_VECTOR(NBIT-1 DOWNTO 0); --Q9.0
			CLK,RESET: IN STD_LOGIC;
			ENABLE: IN STD_LOGIC;
			READY: OUT STD_LOGIC;
			SIGMA: OUT STD_LOGIC_VECTOR(NBIT-2 DOWNTO 0)
		);
END ENTITY WA;

ARCHITECTURE RTL OF WA IS 
SIGNAL Q:UNSIGNED (NBIT-2+INTEGER(CEIL(LOG2(REAL(WINDOW)))) DOWNTO 0);
SIGNAL POSITIVO: UNSIGNED(NBIT-2 DOWNTO 0); --U8.0
SIGNAL TEMP: SIGNED(SAMPLE'RANGE);
SIGNAL MUX: UNSIGNED(NBIT-2 DOWNTO 0);
CONSTANT C1: UNSIGNED(2 DOWNTO 0):=TO_UNSIGNED(5,3); --1.25 U0.2
CONSTANT C2: UNSIGNED(2 DOWNTO 0):=TO_UNSIGNED(6,3); -- 1.58 U0.2 APPROX 6.32
SIGNAL MULCOS:UNSIGNED(C1'RANGE);
SIGNAL AVG_TMP: UNSIGNED (NBIT-2+INTEGER(CEIL(LOG2(REAL(WINDOW))))+C1'LEFT+1 DOWNTO 0); --U17.0
SIGNAL NOISE,NOISE_TMP: UNSIGNED(SIGMA'RANGE);
SIGNAL READY_TMP,READY_TMP_TMP: STD_LOGIC;
SIGNAL CHANGE: STD_LOGIC;
SIGNAL CHANGE_SIGMA : STD_LOGIC;
SIGNAL COUNTER: UNSIGNED(INTEGER(CEIL(LOG2(REAL(WINDOW))))-1 DOWNTO 0);
SIGNAL EN_TMP : STD_LOGIC;



BEGIN


TEMP<=ABS(SIGNED(SAMPLE)); --Q9.0
--POSITIVO<=UNSIGNED(TEMP(NBIT-2 DOWNTO 0)); --U8.0

---IL REGISTRO SERVE PER GARANTIRE UN CORRETTO CONFRONTO TRA GLI INGRESSI E IL SALVATAGGIO IN MEMORIA.
PROCESS(CLK,RESET)
BEGIN
IF(RESET='0') THEN
POSITIVO<=(OTHERS=>'0');
	EN_TMP<='1';
ELSIF(RISING_EDGE(CLK)) THEN
	IF(ENABLE='0') THEN
	POSITIVO<=UNSIGNED(TEMP(NBIT-2 DOWNTO 0));
		EN_TMP<=ENABLE;
	END IF;
END IF;
END PROCESS;



COMPARISON:PROCESS(READY_TMP,POSITIVO,NOISE_TMP)
VARIABLE DIFF: SIGNED(NBIT-1 DOWNTO 0);
BEGIN
	IF(READY_TMP='0') THEN
		DIFF:=SIGNED("0" & POSITIVO)-SIGNED( "0" & NOISE_TMP);
			IF(DIFF(DIFF'HIGH)='1') THEN
				CHANGE<='1';
			ELSE	
				CHANGE<='0';
			END IF;
		ELSE
		 CHANGE<='1';
	END IF;
END PROCESS COMPARISON;
		


MUX<= POSITIVO WHEN (CHANGE='1') ELSE NOISE_TMP;

--SUM<=Q + RESIZE(MUX,SUM'LENGTH);



INDEXING:PROCESS(CLK,RESET)
--VARIABLE COUNTER:UNSIGNED(INTEGER(CEIL(LOG2(REAL(WINDOW))))-1 DOWNTO 0);
BEGIN
IF(RESET='0') THEN
		
		COUNTER<=(OTHERS=>'0');
		READY_TMP<='1';
		CHANGE_SIGMA<='0';
ELSIF (RISING_EDGE(CLK)) THEN
		
	IF(EN_TMP='0') THEN
		
			COUNTER<=COUNTER+1;
			
		IF (COUNTER=WINDOW-1) THEN
			 READY_TMP<='0';
			 COUNTER<=(OTHERS=>'0');
			 CHANGE_SIGMA<='1';
		ELSE
			  CHANGE_SIGMA<='0';
		END IF;
		
	END IF;
	--COUNTERR<=COUNTER;
END IF;
END PROCESS INDEXING;

READY<=READY_TMP;

ACCUMULATIONREG:PROCESS(RESET,CLK)
	BEGIN
		IF(RESET='0') THEN
			Q<=(OTHERS=>'0');
		ELSIF RISING_EDGE(CLK) THEN
				--IF(ENABLE='0') THEN
						IF(CHANGE_SIGMA='1') THEN
							Q<=RESIZE(MUX,Q'LENGTH);
							ELSE
							Q<=Q+RESIZE(MUX,Q'LENGTH);
						END IF;
				--END IF;
		END IF;
END PROCESS ACCUMULATIONREG;


PROCESS(CLK,RESET)
BEGIN
IF(RESET='0') THEN
	READY_TMP_TMP<='1';
ELSIF(RISING_EDGE(CLK)) THEN
	READY_TMP_TMP<=READY_TMP;
END IF;
END PROCESS;



MULCOS<= C1 WHEN (READY_TMP_TMP='1') ELSE C2; 

AVG_TMP<=Q*MULCOS; --U17.0




NOISE_TMP<=AVG_TMP(AVG_TMP'HIGH-1 DOWNTO AVG_TMP'HIGH-9)  WHEN(CHANGE_SIGMA='1') ELSE NOISE;
SIGMA<=STD_LOGIC_VECTOR(NOISE_TMP);

SAVINGSIGMA:PROCESS(CLK,RESET)
	BEGIN
	 IF(RESET='0') THEN
		NOISE<=(OTHERS=>'0');
	ELSIF (RISING_EDGE(CLK)) THEN
			IF(CHANGE_SIGMA='1') THEN
				NOISE<=AVG_TMP(AVG_TMP'HIGH-1 DOWNTO AVG_TMP'HIGH-9);
			END IF;
	END IF;
END PROCESS SAVINGSIGMA;


END ARCHITECTURE RTL;



