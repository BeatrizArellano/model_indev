program shelfmodel_main
  use shelfseas,        only: init_shelfseas, run_shelfseas, end_shelfseas


  implicit none

    call init_shelfseas()   

    call run_shelfseas()

    call end_shelfseas()   
        

end program shelfmodel_main




