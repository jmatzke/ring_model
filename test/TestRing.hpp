#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "SimpleStimulus.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "LuoRudy1991.hpp"

class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    PointStimulus2dCellFactory()
        : AbstractCardiacCellFactory<1>(),
          mpStimulus(new SimpleStimulus(-5e5, 1.0))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        double x = pNode->rGetLocation()[0];

        if (x< 0.2) // ie if x<=0.02 and y<=0.02 (and we are assuming here x,y>=0).
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }

};

class TestSimpleRingRun : public CxxTest::TestSuite
{
public:
    void TestSimpleSimulation()
    {
        HeartConfig::Instance()->SetSimulationDuration(1000.0); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetOutputDirectory("RingRun");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        //HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
        HeartConfig::Instance()->SetVisualizeWithCmgui(true);
        //HeartConfig::Instance()->SetVisualizeWithVtk(true);
        //HeartConfig::Instance()->SetVisualizeWithParallelVtk(true);

        PointStimulus2dCellFactory cell_factory;
/*
	double magnitude = -9.0e3; // uA/cm^2
        double start_time = 50.0;
        double duration = 2; //ms
        HeartConfig::Instance()->SetElectrodeParameters(false, 0, magnitude, start_time, duration);
*/
        BidomainProblem<1> bidomain_problem( &cell_factory );

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.8));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.8));

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400); // 1/cm
        HeartConfig::Instance()->SetCapacitance(1.0); // uF/cm^2

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.01, 1.0);

        bidomain_problem.Initialise();

        bidomain_problem.Solve();

        ReplicatableVector res_repl(bidomain_problem.GetSolution());
        for (unsigned i=0; i<res_repl.GetSize(); i++)
        {
        //    std::cout << res_repl[i] << "\n";
        }

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }
};
