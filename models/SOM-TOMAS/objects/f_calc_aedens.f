
      DOUBLE PRECISION FUNCTION CALC_AEDENS(Mass_ORG,Mass_AMS)

      DOUBLE PRECISION Mass_ORG
      DOUBLE PRECISION Mass_AMS
      DOUBLE PRECISION Mass_TOT
      
      DOUBLE PRECISION rho_org
      PARAMETER(rho_org = 1400.0)

      DOUBLE PRECISION rho_ams
      PARAMETER(rho_ams = 1770.0)

      Mass_TOT = Mass_ORG + Mass_AMS
      
      CALC_AEDENS = Mass_TOT/(Mass_ORG/rho_org + Mass_AMS/rho_ams)

      RETURN
      END FUNCTION
