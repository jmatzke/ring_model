#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "SimpleStimulus.hpp"
#include "RegularStimulusZeroNetCharge.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
#include <cmath>
#define PI 3.14159265

class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    boost::shared_ptr<SimpleStimulus> mpPosStimulus;
    boost::shared_ptr<RegularStimulusZeroNetCharge> mpRegStimulus;
    const int Index;

public:
    PointStimulus2dCellFactory(int InputIndex)
        : AbstractCardiacCellFactory<2>(),
          mpStimulus(new SimpleStimulus(1e4, 15.0)), //blocks birectional conduction initially
          mpPosStimulus(new SimpleStimulus(-1e5, 0.5)), //establishes initial ring current
          mpRegStimulus(new RegularStimulusZeroNetCharge(5e4,8.0,500.0,600.0,900.00)),
          //magnitude, duration, period, start, stop time of biphasic shock
          Index(InputIndex)
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        
        double ThetaOne = (Index + 1)*PI/10.0;
        double XStopOne = 1.65*cos(ThetaOne)+0.01;
        double XStartOne = 1.65*cos(ThetaOne)-0.01;
        double ThetaTwo = ThetaOne + PI;
        double XStopTwo = 1.65*cos(ThetaTwo)+0.01;
        double XStartTwo = 1.65*cos(ThetaTwo)-0.01;
        //std::cout << "XStop: " << XStop << " XStart: " << XStart << std::endl;
        

        if (x<-1.63)
        {
        	if (y>0.0){
            		return new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpStimulus);
            		}
    		else {
    			return new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpPosStimulus);
    			}
        } else if (((x-XStopOne)*(x-XStartOne) <= 0) && y>0) {
        	return new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpRegStimulus);
        }  else if (((x-XStopTwo)*(x-XStartTwo) <= 0) && y<0) {
        	return new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpRegStimulus);
        } else {
            return new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpZeroStimulus);
        }
    }

};

class TestMultiRingRun : public CxxTest::TestSuite
{
private:
	void SingleRun(int RunNumber)
	{
		
		std::stringstream output_dir;
		output_dir << "RingRun/Run_" << RunNumber;
		HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
		HeartConfig::Instance()->SetSimulationDuration(1000.0); //ms
		HeartConfig::Instance()->SetMeshFileName("mesh/test/data/annuli/circular_annulus_480_elements");
		//HeartConfig::Instance()->SetVisualizeWithCmgui(true);
		//HeartConfig::Instance()->SetVisualizeWithVtk(true);
		HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
		HeartConfig::Instance()->SetOutputFilenamePrefix("results");
		HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.5,0.5));
		HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(6.4,6.4));
		HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400); // 1/cm
		HeartConfig::Instance()->SetCapacitance(1.0); // uF/cm^2
		HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.05, 0.1, 1.0);

		PointStimulus2dCellFactory cell_factory(RunNumber);
		BidomainProblem<2> bidomain_problem( &cell_factory );
		bidomain_problem.Initialise();
		bidomain_problem.Solve();
		
	}
public:
    void TestMainSimulation()
    {
        const int NumberOfRuns = 9;
        
	//TODO: Parallelize for-loop
	for (int i=0; i < NumberOfRuns; ++i){
		SingleRun(i);
	}
	/*
        ReplicatableVector res_repl(bidomain_problem.GetSolution());
        for (unsigned i=0; i<res_repl.GetSize(); i++)
        {
            std::cout << res_repl[i] << "\n";
        }
	
        HeartEventHandler::Headings();
        HeartEventHandler::Report();
        */
    }
};
