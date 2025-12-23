program testapp
  use fortuno_serial, only : execute_serial_cmd_app, test_list
  use test_math_angmomentum, only : angmomentum_tests => tests
  use test_math_matrixops, only : matrixops_tests => tests
  implicit none

  call execute_serial_cmd_app(test_list([&
      & matrixops_tests(),&
      & angmomentum_tests()&
    ])&
  )

end program testapp
