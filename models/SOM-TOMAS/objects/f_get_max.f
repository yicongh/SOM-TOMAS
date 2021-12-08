

      DOUBLE PRECISION FUNCTION GET_MAX(SIZE,ARRAY)

      INTEGER          SIZE
      DOUBLE PRECISION ARRAY(SIZE)
      
      INTEGER i
      
      GET_MAX = 0.0
      
      DO i = 1,SIZE
         GET_MAX = MAX(GET_MAX,ABS(ARRAY(i)))
      END DO
      
      END FUNCTION
