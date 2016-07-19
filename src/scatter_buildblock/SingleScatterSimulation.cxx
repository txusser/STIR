 #include "stir/scatter/SingleScatterSimulation.h"

START_NAMESPACE_STIR

const char * const
SingleScatterSimulation::registered_name =
  "Single Scatter Simulation";


SingleScatterSimulation::
SingleScatterSimulation() :
    base_type()
{
    this->set_defaults();
}

SingleScatterSimulation::
~SingleScatterSimulation()
{}


void
SingleScatterSimulation::
initialise_keymap()
{
    this->parser.add_start_key("Single Scatter Simulation Parameters");
    this->parser.add_stop_key("end Single Scatter Simulation Parameters");
}

void
SingleScatterSimulation::
set_defaults()
{

}

void
SingleScatterSimulation::
ask_parameters()
{

}

bool
SingleScatterSimulation::
post_processing()
{

}

std::string
SingleScatterSimulation::
method_info() const
{

}


END_NAMESPACE_STIR

