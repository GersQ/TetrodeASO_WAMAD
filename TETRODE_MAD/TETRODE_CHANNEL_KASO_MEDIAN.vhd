
-- This architecture implements a 4-channel Tetrode arragnment, and ASO with MAD as spike detection algortihm, with a scale factor of 32.


LIBRARY IEEE;
USE IEEE.STD_LOGIC_1164.ALL;
USE IEEE.NUMERIC_STD.ALL;
USE IEEE.MATH_REAL.LOG2;
USE IEEE.MATH_REAL.CEIL;

PACKAGE VHDLTOOL is
TYPE POINTERS IS ARRAY (natural range <>) of UNSIGNED(INTEGER(CEIL(LOG2(REAL(17))))-1 DOWNTO 0);
TYPE PARALLEL IS ARRAY (NATURAL RANGE <>) OF STD_LOGIC_VECTOR(9 DOWNTO 0);
TYPE PARALLEL_IN IS ARRAY (NATURAL RANGE <>) OF STD_LOGIC_VECTOR(10 DOWNTO 0);
TYPE PARALLEL_OUT IS ARRAY (NATURAL RANGE <>) OF STD_LOGIC_VECTOR(9 DOWNTO 0);
TYPE ENERGY_OUT IS ARRAY (NATURAL RANGE <>) OF STD_LOGIC_VECTOR(19 DOWNTO 0);
TYPE SQUARE IS ARRAY (NATURAL RANGE <>) OF STD_LOGIC_VECTOR(22 DOWNTO 0);
TYPE FINALBIT IS ARRAY (NATURAL RANGE <>) OF STD_LOGIC;
end package VHDLTOOL;

LIBRARY IEEE;
USE IEEE.STD_LOGIC_1164.ALL;
USE IEEE.NUMERIC_STD.ALL;
USE IEEE.MATH_REAL.LOG2;
USE IEEE.MATH_REAL.CEIL;
USE WORK.VHDLTOOL.ALL;

ENTITY TETRODE_CHANNEL_KASO_MEDIAN IS 
GENERIC(
			K:INTEGER:=4;
			N:INTEGER:=10;
			PIXEL:INTEGER:=1; --MEAN PIXEL 
			WINDOW: INTEGER:=64;
			ROW_W: INTEGER:=4
			);
PORT(
			SAMPLE	: IN PARALLEL(0 TO ROW_W-1);
			CLK		: IN STD_LOGIC;
			RESET		: IN STD_LOGIC;
			GOOD		: OUT STD_LOGIC;
			AP_BIT	: OUT STD_LOGIC
		);
END ENTITY TETRODE_CHANNEL_KASO_MEDIAN;





ARCHITECTURE RTL OF TETRODE_CHANNEL_KASO_MEDIAN IS 
--MODULES INSTANCES

--------BUTTERWORTH FILTER---------------
COMPONENT BUTTERWORTH_BUFFV IS
GENERIC(
			CHANNEL	: INTEGER; 
			N			: INTEGER 
			); 
			
PORT(
		CLK,RESET			: IN STD_LOGIC;
		READY					: IN STD_LOGIC;
		SAMPLE				: IN STD_LOGIC_VECTOR(N-1 DOWNTO 0);    --QN-1.0
		FILTERED				: OUT STD_LOGIC_VECTOR(N-1 DOWNTO 0); --QN-1.0

		INDEXW,INDEXR_OLD	: IN	UNSIGNED(INTEGER(CEIL(LOG2(REAL(CHANNEL*2))))-1 DOWNTO 0)
		); 
END COMPONENT BUTTERWORTH_BUFFV;

COMPONENT BUTT_CONTROLLER IS
GENERIC (
			CHANNEL: INTEGER
			);
PORT(	
		CLK,RESET, ENABLE	:  IN STD_LOGIC;
		READY					:	OUT STD_LOGIC;
		INDEXW				:  OUT UNSIGNED(INTEGER(CEIL(LOG2(REAL(CHANNEL*2))))-1 DOWNTO 0);
		INDEXR_OLD			:  OUT UNSIGNED(INTEGER(CEIL(LOG2(REAL(CHANNEL*2))))-1 DOWNTO 0)
		);
		
END COMPONENT BUTT_CONTROLLER;

-------KASO--------
COMPONENT KASO_BUFF IS
GENERIC(
			K: INTEGER; 
			N: INTEGER ; 
			PIXEL:INTEGER
			);
PORT(
	  SAMPLE					: IN STD_LOGIC_VECTOR(N-1 DOWNTO 0);
	  READY					: IN STD_LOGIC;
	  CLK,RESET				: IN STD_LOGIC;
	  ENERGY					: OUT STD_LOGIC_VECTOR(2*N-1 DOWNTO 0);
	  INDEXW					: IN UNSIGNED(INTEGER(CEIL(LOG2(REAL((K*PIXEL+1)))))-1 DOWNTO 0);
	  INDEXK					: IN UNSIGNED(INTEGER(CEIL(LOG2(REAL((K*PIXEL+1)))))-1 DOWNTO 0)
	 
	  );
END COMPONENT KASO_BUFF;

COMPONENT KASO_CONTROLLER IS
GENERIC (
			K		: INTEGER; 
			PIXEL	: INTEGER
			);
PORT 	(	
			CLK,RESET,ENABLE	: IN STD_LOGIC;
			READY					: OUT STD_LOGIC;
			INDEXW				: OUT UNSIGNED(INTEGER(CEIL(LOG2(REAL(K*PIXEL+1))))-1 DOWNTO 0);
			INDEXK				: OUT UNSIGNED(INTEGER(CEIL(LOG2(REAL(K*PIXEL+1))))-1 DOWNTO 0)
		);
END COMPONENT KASO_CONTROLLER;

-----SMOOTHING-----------
COMPONENT SMOOTHING_BUFF IS
GENERIC(
			CHANNEL	: INTEGER; 
			K			: INTEGER; 
			N			: INTEGER
			);
PORT( 
		SAMPLE				: IN STD_LOGIC_VECTOR(2*N-1 DOWNTO 0); --Q19.0
		CLK,RESET			: IN STD_LOGIC;
		INDEXR				: IN POINTERS (0 TO 4*K);
		INDEXW				: IN UNSIGNED(INTEGER(CEIL(LOG2(REAL(4*K*CHANNEL+1))))-1 DOWNTO 0);
		READY					: IN STD_LOGIC;
		SMOOTHED				: OUT STD_LOGIC_VECTOR(2*N-1 DOWNTO 0) --Q19.0
		);
END COMPONENT SMOOTHING_BUFF;

COMPONENT SMOOTH_CONTROLLER IS 
GENERIC (
				K:INTEGER; 
				CHANNEL: INTEGER
			);
PORT(
		CLK,RESET,ENABLE	: IN STD_LOGIC;
		READY					: OUT STD_LOGIC;
		INDEXW				: OUT UNSIGNED(INTEGER(CEIL(LOG2(REAL(4*K*CHANNEL+1))))-1 DOWNTO 0);
		INDEXR				: OUT POINTERS (0 TO 4*K)
		);
END COMPONENT SMOOTH_CONTROLLER;



-----------MEDIAN---------------
COMPONENT STATIC_MEDIAN IS 
PORT(
		SAMPLE: IN STD_LOGIC_VECTOR(9 DOWNTO 0);
		READY	: OUT STD_LOGIC;
		CLK,RESET,ENABLE: IN STD_LOGIC;
		MEDIAN: OUT STD_LOGIC_VECTOR(8 DOWNTO 0)
		
		);
END COMPONENT STATIC_MEDIAN;



--SIGNALS DECLARATIONS
SIGNAL MASTER_RESET			: STD_LOGIC;
SIGNAL STARTSAMPLE			: FINALBIT(0 TO ROW_W-1);
SIGNAL STARTINGSAMPLE			: STD_LOGIC;
SIGNAL FILTERED 			: PARALLEL (0 TO ROW_W-1);
SIGNAL INDEXW_BUTT,INDEXR_BUTT		: UNSIGNED(INTEGER(CEIL(LOG2(REAL(PIXEL*2))))-1 DOWNTO 0);
SIGNAL START_MEAN			: STD_LOGIC;
SIGNAL RKNEO				: STD_LOGIC;
SIGNAL INDEXW_KASO,INDEXR_KASO		: UNSIGNED(INTEGER(CEIL(LOG2(REAL(K*PIXEL+1))))-1 DOWNTO 0);
SIGNAL ENERGY				: STD_LOGIC_VECTOR(2*N-1 DOWNTO 0);
SIGNAL INDEXW_S				: UNSIGNED(INTEGER(CEIL(LOG2(REAL(4*K*PIXEL+1))))-1 DOWNTO 0);
SIGNAL INDEXR_S				: POINTERS (0 TO 4*K);
SIGNAL SMOOTHED,SMOOTHED_TMP		: STD_LOGIC_VECTOR(2*N-1 DOWNTO 0);
SIGNAL RSMOOTH				: STD_LOGIC;
SIGNAL START_MEDIAN			: STD_LOGIC;
SIGNAL POS_TEMP				: SIGNED(N-1 DOWNTO 0);
SIGNAL POSITIVO				: STD_LOGIC_VECTOR(N-2 DOWNTO 0);
SIGNAL MEAN_FILTERED_TMP		: SIGNED(N+1 DOWNTO 0);
SIGNAL MEAN_FILTERED			: STD_LOGIC_VECTOR(N-1 DOWNTO 0);
SIGNAL SIGMA				: STD_LOGIC_VECTOR(N-2 DOWNTO 0); 
SIGNAL VARIANCE				: UNSIGNED(2*SIGMA'LENGTH+4 DOWNTO 0); --Q22.0
SIGNAL READY_MEDIAN			: STD_LOGIC;
SIGNAL DIFF				: SIGNED(VARIANCE'RANGE);
SIGNAL CHECK				: STD_LOGIC;

BEGIN

--SYNCHRONIZER RESET
SYNC_RESET:PROCESS(CLK,RESET)
BEGIN
IF(RESET='0') THEN
	MASTER_RESET<='0';

ELSIF(RISING_EDGE(CLK)) THEN
	MASTER_RESET<='1';

END IF;
END PROCESS;

SYNC_start:PROCESS(CLK,RESET)
BEGIN
IF(RESET='0') THEN
	STARTINGSAMPLE<='1';

ELSIF(RISING_EDGE(CLK)) THEN
	if(MASTER_RESET='1') then
	STARTINGSAMPLE<='0';
	end if;
END IF;
END PROCESS;


--BUTTERWORTH FILTER
--BUTTERWORTH FILTER
BUTTERWORTH_CONTROL: BUTT_CONTROLLER 
											GENERIC MAP(
															CHANNEL=>PIXEL
															)
												PORT MAP(
																CLK=>CLK,
																RESET=>MASTER_RESET,
																ENABLE=>STARTINGSAMPLE,
																READY=>START_MEAN,
																INDEXW=>INDEXW_BUTT,
																INDEXR_OLD=>INDEXR_BUTT
															);


	BUTTERWORTH: FOR I IN 0 TO ROW_W-1 GENERATE
					BB:BUTTERWORTH_BUFFV GENERIC MAP(
																	CHANNEL=>PIXEL,
																	N=>N
																)
													PORT MAP(
																	CLK=>CLK,
																	RESET=>MASTER_RESET,
																	SAMPLE=>SAMPLE(I),
																	FILTERED=>FILTERED(I),
																	INDEXW=>INDEXW_BUTT,
																	INDEXR_OLD=>INDEXR_BUTT,
																	READY=>START_MEAN
																	

																);
				  END GENERATE BUTTERWORTH;



	MEAN_FILTERED_TMP<= RESIZE(SIGNED(FILTERED(0)),N+2)
								+RESIZE(SIGNED(FILTERED(1)),N+2)
									+RESIZE(SIGNED(FILTERED(2)),N+2)
										+RESIZE(SIGNED(FILTERED(3)),N+2);
											
	
	MEAN_FILTERED<=STD_LOGIC_VECTOR(MEAN_FILTERED_TMP(N+1 DOWNTO 2)) WHEN (START_MEAN='0') ELSE (OTHERS=>'0');






----KASO							
	KASO_CONTROL: KASO_CONTROLLER
								GENERIC MAP(
													K=>K,
													PIXEL=>PIXEL
												)
									PORT MAP(	
												CLK=>CLK,
												RESET=>MASTER_RESET,
												READY=>RKNEO,
												ENABLE=>START_MEAN,
												INDEXK=>INDEXR_KASO,
												INDEXW=>INDEXW_KASO
												);
			
	


					KASO: KASO_BUFF
							GENERIC MAP(
												K=>K,
												N=>N,
												PIXEL=>PIXEL
											)
								PORT MAP(
												CLK=>CLK,
												RESET=>MASTER_RESET,
												READY=>RKNEO,
												SAMPLE=>MEAN_FILTERED,
												ENERGY=>ENERGY,
												INDEXK=>INDEXR_KASO,
												INDEXW=>INDEXW_KASO
											);
--SMOOTHING	
				SMOOTH:	SMOOTHING_BUFF 
								GENERIC MAP(
													CHANNEL=>PIXEL,
													K=>K,
													N=>N
												)
									PORT MAP(
													CLK=>CLK,
													RESET=>MASTER_RESET,
													INDEXR=>INDEXR_S,
													INDEXW=>INDEXW_S,
													SAMPLE=>ENERGY,
													SMOOTHED=>SMOOTHED_TMP,
													READY=>RSMOOTH
												);

	
	

	
SMOOTH_CONTROL:	   SMOOTH_CONTROLLER						
							GENERIC MAP(
												K=>K,
												CHANNEL=>PIXEL
											)
								PORT MAP(
												CLK=>CLK,
												RESET=>MASTER_RESET,
												ENABLE=>RKNEO,
												INDEXR=>INDEXR_S,
												INDEXW=>INDEXW_S,
												READY=>RSMOOTH
											);
											
----------- SIGMA ESTIMATE	
--ABS
--POS_TEMP<=ABS(SIGNED(MEAN_FILTERED));
--POSITIVO<=STD_LOGIC_VECTOR(POS_TEMP(N-2 DOWNTO 0));


DELAY:PROCESS(CLK,MASTER_RESET)
BEGIN
IF(MASTER_RESET='0') THEN
	START_MEDIAN<='1';
ELSIF(RISING_EDGE(CLK)) THEN
	START_MEDIAN<=START_MEAN;
END IF;
END PROCESS;


MEDIAN_S: STATIC_MEDIAN PORT MAP(
											CLK   =>CLK,
											RESET =>MASTER_RESET,
											ENABLE=>START_MEDIAN,
											READY=>READY_MEDIAN,
											SAMPLE=>MEAN_FILTERED,
											MEDIAN=>SIGMA--MEDIANA HA UN REGISTRO IN USCITA	
											);
	
	
	VARIANCE	<= UNSIGNED(SIGMA)*UNSIGNED(SIGMA) & "00000";
	
	
----COMPARISON SECTION
PROCESS(CLK,MASTER_RESET) 
BEGIN
IF(MASTER_RESET='0') THEN
	SMOOTHED<=(OTHERS=>'0');
ELSIF(RISING_EDGE(CLK)) THEN
	SMOOTHED<=SMOOTHED_TMP;
END IF;
END PROCESS;

GOOD<=READY_MEDIAN;
DIFF<=RESIZE(SIGNED(SMOOTHED),VARIANCE'LENGTH)- RESIZE(SIGNED(VARIANCE),VARIANCE'LENGTH);
CHECK<=DIFF(DIFF'HIGH);
AP_BIT<='1' WHEN (CHECK='0') ELSE '0';


	
	
END ARCHITECTURE RTL;
