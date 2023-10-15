integer*8 function helium_ratio_wrap(argc, argv)
    implicit none

    integer*8 argc, argv(*)

    call helium_ratio_7jun2016(%val(argv(1)), %val(argv(2)), %val(argv(3)), %val(argv(4)), & 
      %val(argv(5)), %val(argv(6)), %val(argv(7)))

    helium_ratio_wrap=1

end
