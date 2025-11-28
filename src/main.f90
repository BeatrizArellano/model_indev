program shelfmodel_main
  use shelfseas,        only: init_shelfseas, run_shelfseas, end_shelfseas
  use sim_timer,        only: sim_timer_start, sim_timer_stop_and_print


  implicit none

    call sim_timer_start()

    call init_shelfseas()   

    call run_shelfseas()

    call end_shelfseas()

    call sim_timer_stop_and_print()       

end program shelfmodel_main




