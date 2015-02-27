#include <Simbody.h>
using namespace SimTK;

namespace {
    const Real g = 9.8;
    const Real m = 1;
    const Real r = 1;
    const Real mr2 = m*r*2;

    class PositionInfo {
    public:
        PositionInfo() : p_AC(NaN), p_BC(NaN), p_FC(NaN) {}
        Vec3 p_AC, p_BC, p_FC; // vectors in ancestor frame
    };
    class VelocityInfo {
    public:
        VelocityInfo() : v_AC(NaN), v_FC(NaN), p_FC_dot(NaN) {}
        Vec3 v_AC, v_FC, p_FC_dot; // vectors in ancestor frame
    };

    class ConstraintRollingDisc: public Constraint::Custom::Implementation {
    public:
        // Constructor takes a plane body, a disc body, and the disc radius.
        // Tell the base class that this constraint generates 1 holonomic
        // (position level), 2 nonholonomic (velocity level), and 0
        // acceleration-only constraint equations.
        ConstraintRollingDisc(
                const SimbodyMatterSubsystem& matter,
                MobilizedBody& planeMobod,
                MobilizedBody& discMobod,
                Real discRadius) :
            Implementation(planeMobod.updMatterSubsystem(), 1, 2, 0),
            m_matterIndex(matter.getMySubsystemIndex()),
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
        virtual void calcPositionErrors
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
        virtual void calcPositionDotErrors
           (const State&                                    state,
            const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB,
            const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
            Array_<Real>&                                   pverr) const
        {
            const Transform& X_AF = getBodyTransformFromState(state, m_plane_F);
            const Transform& X_AB = getBodyTransformFromState(state, m_disc_B);
            const Vec3& w_AF = getBodyAngularVelocity(allV_AB, m_plane_F);
            const Vec3& v_AF = getBodyOriginVelocity(allV_AB, m_plane_F);

            const UnitVec3& Pz_A = X_AF.z(); // m_plane_F normal direction (down) in A

            const Vec3& p_BC_A = getPositionInfo(state).p_BC;

            const Vec3& p_AC = getPositionInfo(state).p_AC;
            const Vec3 v_AC = findStationInAVelocity(state, allV_AB, m_disc_B, p_BC_A);

            const Vec3& p_FC_A = getPositionInfo(state).p_FC;
            const Vec3 p_FC_A_dot = v_AC - v_AF;
            const Vec3 v_FC_A = p_FC_A_dot - w_AF % p_FC_A;

            pverr[0] = ~v_FC_A * Pz_A;
        }

        // Simbody supplies acceleration information in argument list; position and
        // velocity info is in the state. Return second time derivative of position
        // constraint error.
        virtual void calcPositionDotDotErrors
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

            const Vec3& p_BC_A = getPositionInfo(state).p_BC;

            const Vec3& p_AC = getPositionInfo(state).p_AC;
            const Vec3& v_AC = getVelocityInfo(state).v_AC;
            const Vec3 a_AC = findStationInAAcceleration(state, allA_AB, m_disc_B, p_BC_A);

            const Vec3& p_FC_A = getPositionInfo(state).p_FC;
            const Vec3& p_FC_A_dot = getVelocityInfo(state).p_FC_dot;
            const Vec3 p_FC_A_dotdot = a_AC - a_AF;

            const Vec3& v_FC_A = getVelocityInfo(state).v_FC;
            const Vec3 v_FC_A_dot = p_FC_A_dotdot -
                (b_AF % p_FC_A + w_AF % p_FC_A_dot);

            const Vec3 a_FC_A = v_FC_A_dot - w_AF % v_FC_A;

            paerr[0] = ~a_FC_A * Pz_A;
        }

        // Simbody provides calculated constraint multiplier in argument list; we
        // turn that into forces here and apply them to the two bodies at the
        // contact points
        virtual void addInPositionConstraintForces
           (const State&                                state,
            const Array_<Real>&                         multipliers,
            Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForcesInA,
            Array_<Real,      ConstrainedQIndex>&       qForces) const
        {
            const Transform& X_AF = getBodyTransformFromState(state, m_plane_F);
            const UnitVec3& Pz_A = X_AF.z(); // m_plane_F normal direction (down) in A

            const Vec3& p_BC_A = getPositionInfo(state).p_BC;
            const Vec3& p_FC_A = getPositionInfo(state).p_FC;
            const Vec3 force_A = multipliers[0] * Pz_A;

            // manually add in body forces to reduce unnecessary rotations
            bodyForcesInA[m_disc_B] += SpatialVec(p_BC_A % force_A, force_A); // rXf, f
            bodyForcesInA[m_plane_F] -= SpatialVec(p_FC_A % force_A, force_A);
        }


        // Implementation of virtuals required for nonholonomic constraints.
        virtual void calcVelocityErrors
           (const State&                                    state, // Stage::Position
            const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB,
            const Array_<Real,      ConstrainedUIndex>&     constrainedU,
            Array_<Real>&                                   verr) const
        {
            const Transform& X_AF = getBodyTransformFromState(state, m_plane_F);
            const Transform& X_AB = getBodyTransformFromState(state, m_disc_B);
            const Vec3& p_AF = X_AF.p();
            const Vec3& w_AF = getBodyAngularVelocity(allV_AB, m_plane_F);
            const Vec3& v_AF = getBodyOriginVelocity(allV_AB, m_plane_F);

            const Vec3& p_BC_A = getPositionInfo(state).p_BC;

            const Vec3& p_AC = getPositionInfo(state).p_AC;
            const Vec3 v_AC = findStationInAVelocity(state, allV_AB, m_disc_B, p_BC_A);

            const Vec3& p_FC_A = getPositionInfo(state).p_FC;
            const Vec3 p_FC_A_dot = v_AC - v_AF;
            const Vec3 v_FC_A = p_FC_A_dot - w_AF % p_FC_A;

            verr[0] = ~X_AF.x() * v_FC_A; // x-component
            verr[1] = ~X_AF.y() * v_FC_A; // y-component
        }

        virtual void calcVelocityDotErrors
           (const State&                                    state, // Stage::Velocity
            const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB,
            const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
            Array_<Real>&                                   vaerr) const
        {
            const Transform& X_AF = getBodyTransformFromState(state, m_plane_F);
            const Transform& X_AB = getBodyTransformFromState(state, m_disc_B);
            const Vec3& p_AF = X_AF.p();
            const Vec3& w_AF = getBodyAngularVelocityFromState(state, m_plane_F);
            const Vec3& v_AF = getBodyOriginVelocityFromState(state, m_plane_F);
            const Vec3& b_AF = getBodyAngularAcceleration(allA_AB, m_plane_F);
            const Vec3& a_AF = getBodyOriginAcceleration(allA_AB, m_plane_F);

            const Vec3& p_BC_A = getPositionInfo(state).p_BC;

            const Vec3& p_AC = getPositionInfo(state).p_AC;
            const Vec3& v_AC = getVelocityInfo(state).v_AC;
            const Vec3 a_AC = findStationInAAcceleration(state, allA_AB, m_disc_B, p_BC_A);

            const Vec3& p_FC_A = getPositionInfo(state).p_FC;
            const Vec3 p_FC_A_dot = getVelocityInfo(state).p_FC_dot;
            const Vec3 p_FC_A_dotdot = a_AC - a_AF;

            const Vec3 v_FC_A = getVelocityInfo(state).v_FC;
            const Vec3 v_FC_A_dot = p_FC_A_dotdot -
                (b_AF % p_FC_A + w_AF % p_FC_A_dot);

            const Vec3 a_FC_A = v_FC_A_dot - w_AF % v_FC_A;

            vaerr[0] = ~X_AF.x() * a_FC_A; // x-component
            vaerr[1] = ~X_AF.y() * a_FC_A; // y-component
        }

        virtual void addInVelocityConstraintForces
           (const State&                                    state, // Stage::Velocity
            const Array_<Real>&                             multipliers, // mv of these
            Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
            Array_<Real,      ConstrainedUIndex>&           mobilityForces) const
        {
            const Real lambda0 = multipliers[0];
            const Real lambda1 = multipliers[1];

            const Vec3& p_BC_A = getPositionInfo(state).p_BC;
            const Vec3& p_FC_A = getPositionInfo(state).p_FC;

            const Transform& X_AF = getBodyTransformFromState(state, m_plane_F);
            const Vec3 force_A = lambda0*X_AF.x() + lambda1*X_AF.y();

            // manually add in body forces to reduce unnecessary rotations
            bodyForcesInA[m_disc_B] += SpatialVec(p_BC_A % force_A, force_A); // rXf, f
            bodyForcesInA[m_plane_F] -= SpatialVec(p_FC_A % force_A, force_A);
        }

        virtual void realizeTopology(State& state) const {
           auto mThis = const_cast<ConstraintRollingDisc*>(this);
           mThis->m_positionIndex = state.allocateCacheEntry(
               m_matterIndex, Stage::Position, new Value<PositionInfo>());
           mThis->m_velocityIndex = state.allocateCacheEntry(
               m_matterIndex, Stage::Velocity, new Value<VelocityInfo>());
        }

        virtual void realizePosition(const State& state) const {
            realizePositionInfo(state);
        }

        virtual void realizeVelocity(const State& state) const {
            realizeVelocityInfo(state);
        }


    private:
        void realizePositionInfo(const State& state) const {
            if (state.isCacheValueRealized(m_matterIndex, m_positionIndex)) {
                return;
            }
            const Transform& X_AF = getBodyTransformFromState(state, m_plane_F);
            const Transform& X_AB = getBodyTransformFromState(state, m_disc_B);
            const UnitVec3& Pz_A = X_AF.z(); // plane normal direction (down) in A
            const UnitVec3& Dy_A = X_AB.y(); // disc spin axis direction in A
            const UnitVec3 n_A((Dy_A % Pz_A) % Dy_A); // dir from disc origin to contact point

            PositionInfo& pos = updPositionInfo(state);
            pos.p_BC = m_radius * n_A;
            pos.p_AC = pos.p_BC + X_AB.p();
            pos.p_FC = pos.p_AC - X_AF.p();

            state.markCacheValueRealized(m_matterIndex, m_positionIndex);
        }
        const PositionInfo& getPositionInfo(const State& state) const {
            return Value<PositionInfo>::downcast(
                state.getCacheEntry(m_matterIndex, m_positionIndex));
        }
        PositionInfo& updPositionInfo(const State& state) const {
            return Value<PositionInfo>::updDowncast(
                    state.updCacheEntry(m_matterIndex, m_positionIndex));
        }

        void realizeVelocityInfo(const State& state) const {
            if (state.isCacheValueRealized(m_matterIndex, m_velocityIndex)) {
                return;
            }
            const Vec3& w_AF = getBodyAngularVelocityFromState(state, m_plane_F);
            const Vec3& v_AF = getBodyOriginVelocityFromState(state, m_plane_F);
            const Vec3& w_AB = getBodyAngularVelocityFromState(state, m_disc_B);
            const Vec3& v_AB = getBodyOriginVelocityFromState(state, m_disc_B);
            const PositionInfo& pos = getPositionInfo(state);

            const Vec3 p_BC_dot = w_AB % pos.p_BC;

            VelocityInfo& vel = updVelocityInfo(state);
            vel.v_AC = v_AB + p_BC_dot;
            vel.p_FC_dot = vel.v_AC - v_AF;
            vel.v_FC = vel.p_FC_dot - w_AF % pos.p_FC;

            state.markCacheValueRealized(m_matterIndex, m_velocityIndex);
        }
        const VelocityInfo& getVelocityInfo(const State& state) const {
            return Value<VelocityInfo>::downcast(
                state.getCacheEntry(m_matterIndex, m_velocityIndex));
        }
        VelocityInfo& updVelocityInfo(const State& state) const {
            return Value<VelocityInfo>::updDowncast(
                    state.updCacheEntry(m_matterIndex, m_velocityIndex));
        }

        Vec3 findStationInAVelocity
           (const State& state,
            const Array_<SpatialVec, ConstrainedBodyIndex>& allV_AB,
            ConstrainedBodyIndex B,
            const Vec3& p_BS_A) const {
            const Vec3& w_AB = getBodyAngularVelocity(allV_AB, B);
            const Vec3& v_AB = getBodyOriginVelocity(allV_AB, B);

            return v_AB + w_AB % p_BS_A;
        }

        Vec3 findStationInAAcceleration
           (const State& state,
            const Array_<SpatialVec, ConstrainedBodyIndex>& allA_AB,
            ConstrainedBodyIndex B,
            const Vec3& p_BS_A) const {
            const Vec3& w_AB = getBodyAngularVelocityFromState(state, B);
            const Vec3& b_AB = getBodyAngularAcceleration(allA_AB, B);
            const Vec3& a_AB = getBodyOriginAcceleration(allA_AB, B);

            return a_AB + (b_AB % p_BS_A) + w_AB % (w_AB % p_BS_A);
        }

        const SubsystemIndex m_matterIndex;
        ConstrainedBodyIndex m_plane_F;
        ConstrainedBodyIndex m_disc_B;
        const Real m_radius;
        CacheEntryIndex m_positionIndex;
        CacheEntryIndex m_velocityIndex;
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
                    matter, matter.Ground(), disc1, r));

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

        disc1.setDefaultRotation(Rotation(Pi/6, YAxis));
        system.realizeTopology();
        State state = system.getDefaultState();
        disc1.setU(state, Vec6(0, 10, 0, 0, 0, 0));
        //disc2.setU(state, Vec6(-1, 5, 0, 0, 0, 0));

        RungeKuttaMersonIntegrator integ(system);
        TimeStepper ts(system, integ);
        ts.initialize(state);
        ts.stepTo(50.0);
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
