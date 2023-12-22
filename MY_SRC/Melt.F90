!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Mosbeux, Caillet
! *  Email:   cmosbeux@univ-grenoble-alpes.fr
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 31 May 2023
! *  Modification : 24 Novembre 2023 (add of quadratic melting)
! *****************************************************************************/

  
!-------------------------------------------------------------------------------
!> Subroutine for calling different melt parameterizations 
!   * PICO
!   * Quadratic 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE MeltSolver_Init0( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Solver_t), TARGET :: Solver

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='Melt', ParamName

  Params => GetSolverParams()
  ParamName = ListGetString(Params, 'Parameterisation Name', UnFoundFatal=.TRUE.)

  CALL INFO(SolverName,'#------------------------------------------', Level=5)
  CALL INFO(SolverName,'Initialisation Melt Parameterisation.', Level=5)
  CALL INFO(SolverName,'#------------------------------------------', Level=5)

  SELECT CASE (TRIM(ParamName))

    CASE ('pico')
      !------------------------------------------------------------------------
      ! 1 - Check if data has been prescribed in the sif file otherwise set default values
      !------------------------------------------------------------------------
      IF(.NOT. ListCheckPresent(Params,'Bottom Surface Variable Name')) &
        CALL ListAddString(Params,'Bottom Surface Variable Name','zb')
      IF(.NOT. ListCheckPresent(Params,'DistGL Name')) &
        CALL ListAddString(Params,'DistGL Name','distGL')
      IF(.NOT. ListCheckPresent(Params,'DistIF Name')) &
        CALL ListAddString(Params,'DistIF Name','distIF')
      IF(.NOT. ListCheckPresent(Params,'PanAntarctic')) &
        CALL ListAddLogical(Params,'PanAntarctic',.TRUE.)
      IF(.NOT. ListCheckPresent(Params,'Grounding Line Melt')) &
        CALL ListAddLogical(Params,'Grounding Line Melt',.FALSE.)

    CASE ('quadratic')
      !------------------------------------------------------------------------
      ! 1 - Check if data has been prescribed in the sif file otherwise set default values
      !------------------------------------------------------------------------
      IF(.NOT. ListCheckPresent(Params,'X Dim Name')) &
        CALL ListAddString(Params,'X Dim Name','x')
      IF(.NOT. ListCheckPresent(Params,'Y Dim Name')) &
        CALL ListAddString(Params,'Y Dim Name','y')
      IF(.NOT. ListCheckPresent(Params,'Z Dim Name')) &
        CALL ListAddString(Params,'Z Dim Name','z')
      IF(.NOT. ListCheckPresent(Params,'Time Dim Name')) &
        CALL ListAddString(Params,'Time Dim Name','time')
      IF(.NOT. ListCheckPresent(Params,'X Var Name')) &
        CALL ListAddString(Params,'X Var Name','x')
      IF(.NOT. ListCheckPresent(Params,'Y Var Name')) &
        CALL ListAddString(Params,'Y Var Name','y')
      IF(.NOT. ListCheckPresent(Params,'Z Var Name')) &
        CALL ListAddString(Params,'Z Var Name','z')
      IF(.NOT. ListCheckPresent(Params,'Time Var Name')) &
        CALL ListAddString(Params,'Time Var Name','time')
      IF(.NOT. ListCheckPresent(Params,'Elmer Coordinate 3')) &
        CALL ListAddString(Params,'Elmer Coordinate 3','zb')
      IF(.NOT. ListCheckPresent(Params,'Get Cell Value')) &
        CALL ListAddLogical(Params,'Get Cell Value',.TRUE.)
      IF(.NOT. ListCheckPresent(Params,'Allow bound mismatch')) &
        CALL ListAddLogical(Params,'Allow bound mismatch',.TRUE.)
      IF(.NOT. ListCheckPresent(Params,'Read full array')) &
        CALL ListAddLogical(Params,'Read full array',.TRUE.)
      IF(.NOT. ListCheckPresent(Params,'Grounding Line Melt')) &
        CALL ListAddLogical(Params,'Grounding Line Melt',.FALSE.) 
    CASE DEFAULT
        CALL FATAL(SolverName,'Parameterisation does not exist.')
  END SELECT
!------------------------------------------------------------------------------
END SUBROUTINE MeltSolver_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE MeltSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE Distance
  USE PICO
  USE QUADRATIC

  IMPLICIT NONE
  
  TYPE(Model_t)  :: Model
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Solver_t), TARGET :: Solver
  
  REAL(KIND=dp) :: dt, Time, Timestepsize, tpoint, yearinday
  LOGICAL ::  TransientSimulation, UnFoundFatal=.TRUE.
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='Melt', DistGLName, DistIFName, ParamName, FileSalinity, FileTemperature
  INTEGER,SAVE :: VisitedTimes = 0
  INTEGER :: step, Timestep

  Params => GetSolverParams()
  ParamName = ListGetString(Params, 'Parameterisation Name', UnFoundFatal=.TRUE.)
  
  CALL INFO(SolverName,'#------------------------------------', Level=3)
  CALL INFO(SolverName,'Starting Melt Parameterisation.', Level=3)
  CALL INFO(SolverName,'#------------------------------------', Level=3)
  
  VisitedTimes = VisitedTimes + 1

  SELECT CASE (TRIM(ParamName))

    CASE ('pico')
      CALL INFO(Trim(SolverName),'Choice : PICO', Level=3 )
      !------------------------------------
      !1 - Compute distance to GL and IF
      !------------------------------------
      ! Grounding line
      DistGLName = ListGetString(Params, 'DistGL Name', UnFoundFatal=.TRUE.)
      CALL DistanceSolver2( Model,Solver,dt,TransientSimulation,DistGLName)

      ! Ice Front
      CALL IceFrontMask( Model,Solver,dt,TransientSimulation )
      DistIFName = ListGetString(Params, 'DistIF Name', UnFoundFatal=.TRUE.)
      CALL DistanceSolver2( Model,Solver,dt,TransientSimulation,DistIFName)
     
      !---------------------------------------
      !2 - Run the BoxModel Solver (PICO.F90)
      !---------------------------------------
      CALL boxmodel_solver(Model,Solver,dt,TransientSimulation)

    CASE ('quadratic')
      CALL INFO(Trim(SolverName),'Choice : QUADRATIC', Level=3 )
      Time = GetTime()
      Timestep = GetTimestep()
      Timestepsize = GetTimestepsize()
      !PRINT *,'Time verif',Time
      !PRINT *,'NB_Size',Timestep
      !PRINT *,'T_size',Timestepsize
      yearinday = ListGetCReal( Model % Constants, 'Calendar', UnFoundFatal = UnFoundFatal )
      FileSalinity = ListGetString(Params, 'File_Salinity', UnFoundFatal=.TRUE.)
      FileTemperature = ListGetString(Params, 'File_Temperature', UnFoundFatal=.TRUE.)
      !------------------------------------
      !1 - Read 3D file for T and S
      !------------------------------------
      step = NINT(yearinday/Timestepsize) + 1
      IF (VisitedTimes.EQ.1 .OR. MOD(VisitedTimes,step).EQ.0) THEN
      !PRINT *,'step',step
      !PRINT *,'VIsit',VisitedTimes
        tpoint = NINT(Time/yearinday) + 1
        CALL ListAddConstReal(Params,'Index Read',tpoint)
        CALL INFO(SolverName,'#-------------------------------',Level=3)
        CALL INFO(SolverName, 'Index to Read: '//I2S(NINT(tpoint)),Level=3)
        CALL INFO(SolverName,'#-------------------------------',Level=3)
        
        ! Salinity
        CALL ListAddString(Params,'Filename',FileSalinity,CaseConversion=.FALSE.)
        CALL ListAddString(Params,'Variable 1','salinity')
        CALL ListAddString(Params,'Target Variable 1','salinity')
        CALL GridDataReader( Model,Solver,dt,TransientSimulation )
        ! Temperature
        CALL ListAddString(Params,'Filename',FileTemperature,CaseConversion=.FALSE.)
        CALL ListAddString(Params,'Variable 1','temperature')
        CALL ListAddString(Params,'Target Variable 1','temperature')
        CALL GridDataReader( Model,Solver,dt,TransientSimulation )
      END IF
      !---------------------------------------
      !2 - Run the quadratic Solver (QUADRATIC.F90)
      !---------------------------------------
      CALL quadratic_solver( Model,Solver,dt,TransientSimulation )
      
    CASE DEFAULT
        CALL FATAL(SolverName,'Parameterisation does not exist.')
  END SELECT  

END SUBROUTINE MeltSolver
