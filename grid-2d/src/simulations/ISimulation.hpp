#pragma once

namespace fsim {
    class ISimulation {
        public:
        virtual void step() = 0;
    };
}
