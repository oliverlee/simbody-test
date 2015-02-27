#include <Simbody.h>
using namespace SimTK;

namespace {
    const Real g = 9.8;
    const Real m = 1;
    const Real r = 1;
    const Real mr2 = m*r*2;

    class ConstraintRollingDisc: public Constraint::Custom::Implementation {
    public:
        // Constructor takes a plane body, a disc body, and the disc radius.
        // Tell the base class that this constraint generates 1 holonomic
        // (position level), 2 nonholonomic (velocity level), and 0
        // acceleration-only constraint equations.
        ConstraintRollingDisc(
                MobilizedBody& planeMobod,
                MobilizedBody& discMobod,
                Real discRadius) :
            Implementation(planeMobod.updMatterSubsystem(), 1, 0, 0),
            m_radius(discRadius) {
                m_plane_F = addConstrainedBody(planeMobod);
                m_disc_B = addConstrainedBody(discMobod);
        }

        // Implement required pure virtual method.
        Implementation* clone () const {return new ConstraintRollingDisc(*this);}

        // Implement the Implementation virtuals required for a holonomic
        // (position level) constraint.

        // Simbody supplies position information in argument list; we calculate
        // the constraint error that represents here.
        void calcPositionErrors
           (const State&                                    state,
            const Array_<Transform,ConstrainedBodyIndex>&   allX_AB,
            const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
            Array_<Real>&                                   perr) const
        {
            const Transform& X_AF = getBodyTransform(allX_AB, m_plane_F);
            const Transform& X_AB = getBodyTransform(allX_AB, m_disc_B);

            const UnitVec3& Pz_A = X_AF.z(); // m_plane_F normal direction (down) in A
            const UnitVec3& Dy_A = X_AB.y(); // m_disc_B direction in A
            const UnitVec3 n_A((Dy_A % Pz_A) % Dy_A); // vector pointing from the m_disc_B origin to m_plane_F contact
            const Vec3 p_BC_A = m_radius * n_A;

            // location of m_disc_B contact point relative to m_plane_F origin
            const Vec3 p_FC_A = X_AB.p() + p_BC_A - X_AF.p();

            // look at only the z component in the m_plane_F frame
            perr[0] = ~p_FC_A * Pz_A;
        }

        // Simbody supplies velocity information in argument list; position info
        // is in the state. Return time derivative of position constraint error.
        void calcPositionDotErrors
           (const State&                                    state,
            const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB,
            const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
            Array_<Real>&                                   pverr) const
        {
            const Transform& X_AF = getBodyTransformFromState(state, m_plane_F);
            const Transform& X_AB = getBodyTransformFromState(state, m_disc_B);
            const Vec3& p_AF = X_AF.p();
            const Vec3& w_AF = getBodyAngularVelocity(allV_AB, m_plane_F);
            const Vec3& v_AF = getBodyOriginVelocity(allV_AB, m_plane_F);

            const UnitVec3& Pz_A = X_AF.z(); // m_plane_F normal direction (down) in A
            const UnitVec3& Dy_A = X_AB.y(); // m_disc_B direction in A
            const UnitVec3 n_A((Dy_A % Pz_A) % Dy_A);
            const Vec3 p_BC_B = X_AB.RInv() * m_radius * n_A;

            const Vec3 p_AO = findStationLocationFromState(state, m_disc_B, p_BC_B);
            const Vec3 v_AO = findStationVelocity(state, allV_AB, m_disc_B, p_BC_B);

            const Vec3 p_FC_A = p_AO - p_AF;
            const Vec3 p_FC_A_dot = v_AO - v_AF;
            const Vec3 v_FC_A = p_FC_A_dot - w_AF % p_FC_A;

            pverr[0] = ~v_FC_A * Pz_A;
        }

        // Simbody supplies acceleration information in argument list; position and
        // velocity info is in the state. Return second time derivative of position
        // constraint error.
        void calcPositionDotDotErrors
           (const State&                                    state,
            const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB,
            const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
            Array_<Real>&                                   paerr) const
        {
            const Transform& X_AF = getBodyTransformFromState(state, m_plane_F);
            const Transform& X_AB = getBodyTransformFromState(state, m_disc_B);
            const Vec3& p_AF = X_AF.p();
            const Vec3& w_AF = getBodyAngularVelocityFromState(state, m_plane_F);
            const Vec3& v_AF = getBodyOriginVelocityFromState(state, m_plane_F);
            const Vec3& b_AF = getBodyAngularAcceleration(allA_AB, m_plane_F);
            const Vec3& a_AF = getBodyOriginAcceleration(allA_AB, m_plane_F);

            const UnitVec3& Pz_A = X_AF.z(); // m_plane_F normal direction (down) in A
            const UnitVec3& Dy_A = X_AB.y(); // m_disc_B direction in A
            const UnitVec3 n_A((Dy_A % Pz_A) % Dy_A);
            const Vec3 p_BC_B = X_AB.RInv() * m_radius * n_A;

            const Vec3 p_AO = findStationLocationFromState(state, m_disc_B, p_BC_B);
            const Vec3 v_AO = findStationVelocityFromState(state, m_disc_B, p_BC_B);
            const Vec3 a_AO = findStationAcceleration(state, allA_AB, m_disc_B, p_BC_B);

            const Vec3 p_FC_A = p_AO - p_AF;
            const Vec3 p_FC_A_dot = v_AO - v_AF;
            const Vec3 p_FC_A_dotdot = a_AO - a_AF;

            const Vec3 v_FC_A = p_FC_A_dot - w_AF % p_FC_A;
            const Vec3 v_FC_A_dot = p_FC_A_dotdot -
                (b_AF % p_FC_A + w_AF % p_FC_A_dot);

            const Vec3 a_FC_A = v_FC_A_dot - w_AF % v_FC_A;

            paerr[0] = ~a_FC_A * Pz_A;
        }

        // Simbody provides calculated constraint multiplier in argument list; we
        // turn that into forces here and apply them to the two bodies at the
        // contact points
        void addInPositionConstraintForces
           (const State&                                state,
            const Array_<Real>&                         multipliers,
            Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForcesInA,
            Array_<Real,      ConstrainedQIndex>&       qForces) const
        {
            const Transform& X_AF = getBodyTransformFromState(state, m_plane_F);
            const Transform& X_AB = getBodyTransformFromState(state, m_disc_B);
            const UnitVec3& Pz_A = X_AF.z(); // m_plane_F normal direction (down) in A
            const UnitVec3& Dy_A = X_AB.y(); // m_disc_B direction in A
            const UnitVec3 n_A((Dy_A % Pz_A) % Dy_A);

            const Vec3 p_BC_A = m_radius * n_A;
            const Vec3 p_FC_A = p_BC_A + X_AB.p() - X_AF.p();

            const Vec3 p_BC_B = X_AB.RInv() * p_BC_A;
            const Vec3 p_FC_F = X_AF.RInv() * p_FC_A;
            //p_FC_F[ZAxis] = 0;

            const Vec3 force_A = multipliers[0] * Pz_A;
            addInStationForce(state, m_disc_B, p_BC_B,  force_A, bodyForcesInA);
            addInStationForce(state, m_plane_F, p_FC_F, -force_A, bodyForcesInA);
        }

    private:
        ConstrainedBodyIndex m_plane_F;
        ConstrainedBodyIndex m_disc_B;
        Real m_radius;
    };
} // namespace

int main() {
    try {
        MultibodySystem system;
        system.setUpDirection(-ZAxis);
        SimbodyMatterSubsystem matter(system);
        GeneralForceSubsystem forces(system);
        Force::UniformGravity gravity(forces, matter, Vec3(0, 0, g));

        Body::Rigid discBody(MassProperties(
            m, Vec3(0), Inertia(mr2/4, mr2/4, mr2/4)));
        // Rotate displayed circle so that is lies in the x-z plane instead of
        // the x-y plane.
        discBody.addDecoration(Rotation(Pi/2, XAxis), DecorativeCircle(r));

        MobilizedBody::Free disc1(matter.Ground(),
                Transform(), discBody, Vec3(0, 0, r));
        Constraint::Custom roll1(new ConstraintRollingDisc(
                    matter.Ground(), disc1, r));

        //MobilizedBody::Free disc2(matter.Ground(),
        //    Transform(), discBody, Vec3(2, 0, r));
        //Constraint::Custom roll2(new ConstraintRollingDisc(
        //            matter.Ground(), disc2, r));

        Visualizer viz(system);
        viz.setShowFrameRate(true);
        viz.setShowSimTime(true);
        viz.setShowFrameNumber(true);
        viz.addFrameController(new Visualizer::BodyFollower(disc1, Vec3(0),
                    Vec3(1, 1, -1)*5, -ZAxis));

        system.addEventReporter(new Visualizer::Reporter(viz, 0.01));

        disc1.setDefaultRotation(Rotation(Pi/4, YAxis));
        system.realizeTopology();
        State state = system.getDefaultState();
        disc1.setU(state, Vec6(1, 0, 5, 0, 0, 0));
        //disc2.setU(state, Vec6(-1, 5, 0, 0, 0, 0));

        RungeKuttaMersonIntegrator integ(system);
        TimeStepper ts(system, integ);
        ts.initialize(state);
        ts.stepTo(10.0);
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
