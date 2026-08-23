program testapp
  use fortuno_serial, only : execute_serial_cmd_app, test_list
  use test_geometry_neighbours, only : neighbours_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      neighbours_tests()&
    ])&
  )

end program testapp
