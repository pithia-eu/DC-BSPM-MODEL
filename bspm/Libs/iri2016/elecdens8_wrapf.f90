integer*8 function elecdens8_wrapf(argc, argv)
    implicit none

    integer*8 argc, argv(*)

    call elecdens7jun2016(%val(argv(1)), %val(argv(2)), %val(argv(3)), %val(argv(4)), &
      %val(argv(5)), %val(argv(6)), %val(argv(7)))

    elecdens8_wrapf=1

end
