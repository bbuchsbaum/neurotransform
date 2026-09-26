test_that("spherical plans match independent Workbench BARYCENTRIC weights", {
  skip_if_not_installed("jsonlite")
  fixture <- jsonlite::read_json(system.file("extdata", "barycentric_oracle", "oracle.json",
                                            package = "neurotransform"), simplifyVector = TRUE)
  expect_identical(fixture$method, "BARYCENTRIC")
  for (i in seq_len(nrow(fixture$cases))) {
    case <- fixture$cases[i, ]
    vertices <- case$vertices[[1]]
    query <- case$query[[1]]
    faces <- case$faces[[1]] + 1L
    expected <- case$weights[[1]]
    for (order in list(seq_len(nrow(faces)), rev(seq_len(nrow(faces))))) {
      plan <- surface_resampling_plan(surface_mesh(query), surface_mesh(vertices, faces[order, ]))
      # No application normalization to hide missing mass in compiled triplets.
      got <- apply_surface_resampling(plan, diag(nrow(vertices)), normalize = "none")
      expect_true(max(abs(got - expected)) < fixture$tolerance, info = case$name)
      expect_equal(rowSums(got), rep(1, nrow(query)), tolerance = 1e-12)
    }
  }
})
