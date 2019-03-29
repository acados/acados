<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile>
  <compound kind="file">
    <name>dense_qp_common.h</name>
    <path>/home/tom/sources/syscop/acados/acados/dense_qp/</path>
    <filename>d8/d1f/dense__qp__common_8h</filename>
    <class kind="struct">qp_solver_config</class>
    <class kind="struct">dense_qp_info</class>
    <member kind="define">
      <type>#define</type>
      <name>QP_SOLVER_CONFIG_</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a1fa43f93fb388164758201bf73ecc8d3</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct d_dense_qp_dim</type>
      <name>dense_qp_dims</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a1fa3f9aada3b9689d3c056c4250a4e47</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct d_dense_qp</type>
      <name>dense_qp_in</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a2929602d12d83114bd9f913d05b7b512</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct d_dense_qp_sol</type>
      <name>dense_qp_out</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a240e134b8710f3880712258654c78014</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct d_dense_qp_res</type>
      <name>dense_qp_res</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a901f48eb561581b958c557c66114ce58</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct d_dense_qp_res_workspace</type>
      <name>dense_qp_res_ws</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>af385c0adced3afb164c620589e8c1004</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>dense_qp_solver_config_calculate_size</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a518de8e97dac89afb208724edc2802ea</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>qp_solver_config *</type>
      <name>dense_qp_solver_config_assign</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a89d6e161183159ed23452e80105211de</anchor>
      <arglist>(void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>dense_qp_dims_calculate_size</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a1559fb7c9ea0a9337b08df8aef79a4ae</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>dense_qp_dims *</type>
      <name>dense_qp_dims_assign</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>aad8e975be4dd617a0e277d495961d7a3</anchor>
      <arglist>(void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>dense_qp_dims_set</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a8405f7f45252221a7741708160871d59</anchor>
      <arglist>(void *config_, void *dims_, const char *field, const int *value)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>dense_qp_in_calculate_size</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a9e6a2b4d85602976ef3c2e208ad5fb95</anchor>
      <arglist>(void *config, dense_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>dense_qp_in *</type>
      <name>dense_qp_in_assign</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a76e03e96b85175253b84362ebf0d0874</anchor>
      <arglist>(void *config, dense_qp_dims *dims, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>dense_qp_out_calculate_size</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>ab4b6935840173d2a5150b91b905fe76a</anchor>
      <arglist>(void *config, dense_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>dense_qp_out *</type>
      <name>dense_qp_out_assign</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a77e741f7374b7c8b057872bff7aed6b3</anchor>
      <arglist>(void *config, dense_qp_dims *dims, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>dense_qp_res_calculate_size</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>ac624bbde45761e42759867e08c57a2a0</anchor>
      <arglist>(dense_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>dense_qp_res *</type>
      <name>dense_qp_res_assign</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a7248b9ac2390098f68fe77e99f74d82a</anchor>
      <arglist>(dense_qp_dims *dims, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>dense_qp_res_workspace_calculate_size</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a47797a1980070301fa888223a542a88b</anchor>
      <arglist>(dense_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>dense_qp_res_ws *</type>
      <name>dense_qp_res_workspace_assign</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a3ba3643f502b3b325bfe39d47ffa5075</anchor>
      <arglist>(dense_qp_dims *dims, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>dense_qp_compute_t</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>abb029e8ffa1ffc87ebb3245de6c7138e</anchor>
      <arglist>(dense_qp_in *qp_in, dense_qp_out *qp_out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>dense_qp_res_compute</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a59006f94a1b9bd5d16891da7f7244ace</anchor>
      <arglist>(dense_qp_in *qp_in, dense_qp_out *qp_out, dense_qp_res *qp_res, dense_qp_res_ws *res_ws)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>dense_qp_res_compute_nrm_inf</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a0ff8e4b9f0f95f86b7849c2509d14a0d</anchor>
      <arglist>(dense_qp_res *qp_res, double res[4])</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>dense_qp_stack_slacks_dims</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>ae75d12fe8225423707f2743eb4b03030</anchor>
      <arglist>(dense_qp_dims *in, dense_qp_dims *out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>dense_qp_stack_slacks</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>aaa3a2ea565150fda2a8b2159de30767c</anchor>
      <arglist>(dense_qp_in *in, dense_qp_in *out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>dense_qp_unstack_slacks</name>
      <anchorfile>d8/d1f/dense__qp__common_8h.html</anchorfile>
      <anchor>a9a256cb1bc2eec8346f082e2aeccf44b</anchor>
      <arglist>(dense_qp_out *in, dense_qp_in *qp_out, dense_qp_out *out)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>ocp_qp_common.h</name>
    <path>/home/tom/sources/syscop/acados/acados/ocp_qp/</path>
    <filename>d5/d0d/ocp__qp__common_8h</filename>
    <class kind="struct">qp_solver_config</class>
    <class kind="struct">ocp_qp_condensing_config</class>
    <class kind="struct">ocp_qp_xcond_solver_config</class>
    <class kind="struct">ocp_qp_info</class>
    <member kind="define">
      <type>#define</type>
      <name>QP_SOLVER_CONFIG_</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a1fa43f93fb388164758201bf73ecc8d3</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct d_ocp_qp_dim</type>
      <name>ocp_qp_dims</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a7ef2c9384130185c57bb2131efb101d9</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct d_ocp_qp</type>
      <name>ocp_qp_in</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>ae4016531493126f5de5a2246087da082</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct d_ocp_qp_sol</type>
      <name>ocp_qp_out</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a7af9d9e132438d87d3eb8b46bac67599</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct d_ocp_qp_res</type>
      <name>ocp_qp_res</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a6fabc999232c68d04d3f6cd1d563df8b</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>struct d_ocp_qp_res_workspace</type>
      <name>ocp_qp_res_ws</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a0a02cefb562d3fcea0656b9bb935a91d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_solver_config_calculate_size</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a92ee0f452537db760cefa9555074cdf8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>qp_solver_config *</type>
      <name>ocp_qp_solver_config_assign</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>aeaa4d97b0fd75ad4928adc261268e2ab</anchor>
      <arglist>(void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_xcond_solver_config_calculate_size</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>ad03804b3cda8d8ea5397ecb2d8a2728a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_xcond_solver_config *</type>
      <name>ocp_qp_xcond_solver_config_assign</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a2f7a271b45b07844467e57c82f8d7c20</anchor>
      <arglist>(void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_condensing_config_calculate_size</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a829662c5823aa1f204e4d4d6b3a22f4e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_condensing_config *</type>
      <name>ocp_qp_condensing_config_assign</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>abf506b8212d5d4e5f759a606472a66b0</anchor>
      <arglist>(void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_dims_calculate_size</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>aef525d20a93aab5d28890d4cf645c295</anchor>
      <arglist>(int N)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_dims *</type>
      <name>ocp_qp_dims_assign</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>ad48b7ba4d4623d5df1ae8f76427d78c9</anchor>
      <arglist>(int N, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_dims_set</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a99bf0e7d9ac23460770b68b2b03149c2</anchor>
      <arglist>(void *config_, void *dims_, int stage, const char *field, const int *value)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_in_calculate_size</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>ae80942aa51701b346c6d19de7ae07ac2</anchor>
      <arglist>(void *config, ocp_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_in *</type>
      <name>ocp_qp_in_assign</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a84ea4f66a7bf3fefe1d46974fed17819</anchor>
      <arglist>(void *config, ocp_qp_dims *dims, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_out_calculate_size</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a2392b33aec29746325c32eb5a549e272</anchor>
      <arglist>(void *config, ocp_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_out *</type>
      <name>ocp_qp_out_assign</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a8369ddbe912d7f13cca3660c0430099d</anchor>
      <arglist>(void *config, ocp_qp_dims *dims, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_res_calculate_size</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a7dfd519601b969cbf089de219e89692d</anchor>
      <arglist>(ocp_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_res *</type>
      <name>ocp_qp_res_assign</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a757b09654e19bc4ac1c6e4e076f0f85f</anchor>
      <arglist>(ocp_qp_dims *dims, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_res_workspace_calculate_size</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a32e4617e0b919d894f99657d927a710c</anchor>
      <arglist>(ocp_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_res_ws *</type>
      <name>ocp_qp_res_workspace_assign</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>af5ba37c389d29bf2353668d6f9bc12ce</anchor>
      <arglist>(ocp_qp_dims *dims, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_res_compute</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a020884f9077247fa7c8ff22221b7c360</anchor>
      <arglist>(ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_qp_res *qp_res, ocp_qp_res_ws *res_ws)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_res_compute_nrm_inf</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a0b91ce68864328f4c49e6ef521331da4</anchor>
      <arglist>(ocp_qp_res *qp_res, double res[4])</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_stack_slacks_dims</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>af4b7f2d906a3cf6fd78dc3dbbbbe18a1</anchor>
      <arglist>(ocp_qp_dims *in, ocp_qp_dims *out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_stack_slacks</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a5f0cabd84948396cbce25940bc8bec9a</anchor>
      <arglist>(ocp_qp_in *in, ocp_qp_in *out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_compute_t</name>
      <anchorfile>d5/d0d/ocp__qp__common_8h.html</anchorfile>
      <anchor>a955481238743f1bea845134c3a65400e</anchor>
      <arglist>(ocp_qp_in *qp_in, ocp_qp_out *qp_out)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>sim_common.h</name>
    <path>/home/tom/sources/syscop/acados/acados/sim/</path>
    <filename>dc/d3b/sim__common_8h</filename>
    <includes id="dd/ddf/external__function__generic_8h" name="external_function_generic.h" local="yes" imported="no">acados/utils/external_function_generic.h</includes>
    <class kind="struct">sim_in</class>
    <class kind="struct">sim_info</class>
    <class kind="struct">sim_out</class>
    <class kind="struct">sim_opts</class>
    <class kind="struct">sim_config</class>
    <member kind="define">
      <type>#define</type>
      <name>NS_MAX</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>a5b4aee33d4791c4dda5d2b1e947db031</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumeration">
      <type></type>
      <name>sim_function_t</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>ad80863c1322e7ac0a828a4273ebdf86d</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>EXPL_ODE_FUN</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>ad80863c1322e7ac0a828a4273ebdf86da7e92ea3e6cf66536ce83bcc12924a029</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>EXPL_ODE_HES</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>ad80863c1322e7ac0a828a4273ebdf86da352bfaa8be199397ab4ee93e4b504c2e</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>EXPL_VDE_FOR</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>ad80863c1322e7ac0a828a4273ebdf86da1201ba5ab145f57451373a86c1ea12fb</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>EXPL_VDE_ADJ</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>ad80863c1322e7ac0a828a4273ebdf86da72a886281193f280a812e18a7fa67505</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>IMPL_ODE_FUN</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>ad80863c1322e7ac0a828a4273ebdf86da3b72e36db7422d0b36d8ad94935ae995</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>IMPL_ODE_FUN_JAC_X_XDOT</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>ad80863c1322e7ac0a828a4273ebdf86da803af88deae80a2ae3e4f0f1dd22731b</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>IMPL_ODE_JAC_X_XDOT_U</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>ad80863c1322e7ac0a828a4273ebdf86da264c7393d37e58ae7f9c24d581b1364c</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>IMPL_ODE_FUN_JAC_X_XDOT_U</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>ad80863c1322e7ac0a828a4273ebdf86da68b5bb721357d1f583e1b6c669cb01df</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>IMPL_ODE_HESS</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>ad80863c1322e7ac0a828a4273ebdf86da786f4a538d7f63bfa3b83bac899aca4f</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>PHI_FUN</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>ad80863c1322e7ac0a828a4273ebdf86da826e8fd394a78af775f3488ee59f2643</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>PHI_FUN_JAC_Y</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>ad80863c1322e7ac0a828a4273ebdf86dac0249e59f9be6a9bcf4773b0672d76af</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>PHI_JAC_Y_UHAT</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>ad80863c1322e7ac0a828a4273ebdf86da5e1945d0baa8f2b9e5ea56d6b36524eb</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>LO_FUN</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>ad80863c1322e7ac0a828a4273ebdf86da82a4e62a8cc8580f6b25a8ffff5e6945</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>GET_GNSF_MATRICES</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>ad80863c1322e7ac0a828a4273ebdf86daa98d5bc721dc1f3d96524aee793ffa87</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_config_calculate_size</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>affd708a09a6abdc4e415bdda0109fa70</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>sim_config *</type>
      <name>sim_config_assign</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>aa4849b52e62fb9f5adb79d8ba0ea13a5</anchor>
      <arglist>(void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_in_calculate_size</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>a6f310f5d179d8dc91c06ed8e7555864d</anchor>
      <arglist>(void *config, void *dims)</arglist>
    </member>
    <member kind="function">
      <type>sim_in *</type>
      <name>sim_in_assign</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>a28a02a3c2b18fd917e497dfa0d9a1acd</anchor>
      <arglist>(void *config, void *dims, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_in_set_</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>a2e6dc4de23ed73a25f525e8d1c1ebf85</anchor>
      <arglist>(void *config_, void *dims_, sim_in *in, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_out_calculate_size</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>a9477c7543e5da99085787a9cf459a739</anchor>
      <arglist>(void *config, void *dims)</arglist>
    </member>
    <member kind="function">
      <type>sim_out *</type>
      <name>sim_out_assign</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>a472b451a2c109938c5796b2d8ca7f8b5</anchor>
      <arglist>(void *config, void *dims, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_out_get_</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>a4fd9c5064c97185131691bf531164a90</anchor>
      <arglist>(void *config, void *dims, sim_out *out, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_opts_set_</name>
      <anchorfile>dc/d3b/sim__common_8h.html</anchorfile>
      <anchor>a453df9f416db1c852646cc15a5226e95</anchor>
      <arglist>(sim_opts *opts, const char *field, void *value)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>external_function_generic.h</name>
    <path>/home/tom/sources/syscop/acados/acados/utils/</path>
    <filename>dd/ddf/external__function__generic_8h</filename>
    <class kind="struct">colmaj_args</class>
    <class kind="struct">blasfeo_dmat_args</class>
    <class kind="struct">blasfeo_dvec_args</class>
    <class kind="struct">external_function_generic</class>
    <class kind="struct">external_function_casadi</class>
    <class kind="struct">external_function_param_casadi</class>
    <member kind="enumeration">
      <type></type>
      <name>ext_fun_arg_t</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a0e56107b5045e5883897e3c5ee7574bc</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>COLMAJ</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a0e56107b5045e5883897e3c5ee7574bca87fe749755e3f9b36d6ad00bacd09394</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>BLASFEO_DMAT</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a0e56107b5045e5883897e3c5ee7574bca522af3b7e0644e3aeaeca02d022759db</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>BLASFEO_DVEC</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a0e56107b5045e5883897e3c5ee7574bca8d7c0b41e50dc79649e3155dfcf53a8b</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>COLMAJ_ARGS</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a0e56107b5045e5883897e3c5ee7574bca1a8585f7653d1edada077b757c674c39</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>BLASFEO_DMAT_ARGS</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a0e56107b5045e5883897e3c5ee7574bca52ef324f1cc736e2bba92f216fbdc82d</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>BLASFEO_DVEC_ARGS</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a0e56107b5045e5883897e3c5ee7574bcad0b4f41c98ede492878c051ba21af873</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>IGNORE_ARGUMENT</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a0e56107b5045e5883897e3c5ee7574bca09823fd544e199c8cd92fe6eed4d3c83</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>external_function_casadi_struct_size</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a67f214b9a3819fb09fd2c10876d75ff2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_set_fun</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>ab1ef4de7baf026fb52990c9728bb6ff1</anchor>
      <arglist>(external_function_casadi *fun, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_set_work</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a05e41c0b365c20851b4c783c39c1a257</anchor>
      <arglist>(external_function_casadi *fun, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_set_sparsity_in</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a90eea79eff1ceef63c7ab895e1b1eee8</anchor>
      <arglist>(external_function_casadi *fun, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_set_sparsity_out</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a49801ce99149dfd375e53067054f2999</anchor>
      <arglist>(external_function_casadi *fun, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_set_n_in</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>ac29321a2fcca3dcd9cbcdc36243a0417</anchor>
      <arglist>(external_function_casadi *fun, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_set_n_out</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>abb1e8d7a9c6880e21b7a69450673e27c</anchor>
      <arglist>(external_function_casadi *fun, void *value)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>external_function_casadi_calculate_size</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a7f8adfb56c89ca394442937cbfc31755</anchor>
      <arglist>(external_function_casadi *fun)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_assign</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>abe97a57bbb51d53bc442a7d575ec236b</anchor>
      <arglist>(external_function_casadi *fun, void *mem)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_wrapper</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>ab8df131720b7ee4775145ff6650fdf61</anchor>
      <arglist>(void *self, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>external_function_param_casadi_struct_size</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>ac711606b6ce47da0b83208463aa37c6a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_set_fun</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a368b879bb5453eb744ddef6cbee17dbe</anchor>
      <arglist>(external_function_param_casadi *fun, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_set_work</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a90ff23432aa5837875727ddaf1dd3f4d</anchor>
      <arglist>(external_function_param_casadi *fun, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_set_sparsity_in</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a5a94d9c87091f3d61b844d5fd4d07919</anchor>
      <arglist>(external_function_param_casadi *fun, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_set_sparsity_out</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a265bcd8e751f8a1ce4c02a223a6990d8</anchor>
      <arglist>(external_function_param_casadi *fun, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_set_n_in</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>aa8daddd2c0e2256c6dfd98256c98be63</anchor>
      <arglist>(external_function_param_casadi *fun, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_set_n_out</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a18b9b218701175392cb3de50e451189b</anchor>
      <arglist>(external_function_param_casadi *fun, void *value)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>external_function_param_casadi_calculate_size</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>afb6ce1af2de59aef7169139dfdebf0ae</anchor>
      <arglist>(external_function_param_casadi *fun, int np)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_assign</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a4bca3685b4070277e5d7471b081a1169</anchor>
      <arglist>(external_function_param_casadi *fun, void *mem)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_wrapper</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>a4f878080d085e516ab3a0b646d4a6244</anchor>
      <arglist>(void *self, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_set_param</name>
      <anchorfile>dd/ddf/external__function__generic_8h.html</anchorfile>
      <anchor>af0743e801143e306eaea0b5934e6fb3f</anchor>
      <arglist>(void *self, double *p)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>condensing_interface.c</name>
    <path>/home/tom/sources/syscop/acados/interfaces/acados_c/</path>
    <filename>d6/dde/condensing__interface_8c</filename>
    <includes id="df/d4f/condensing__interface_8h" name="condensing_interface.h" local="yes" imported="no">acados_c/condensing_interface.h</includes>
    <member kind="function">
      <type>ocp_qp_condensing_config *</type>
      <name>ocp_qp_condensing_config_create</name>
      <anchorfile>d6/dde/condensing__interface_8c.html</anchorfile>
      <anchor>ae3a310a44fe882a61de503d27163c4c7</anchor>
      <arglist>(condensing_plan *plan)</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>ocp_qp_condensing_opts_create</name>
      <anchorfile>d6/dde/condensing__interface_8c.html</anchorfile>
      <anchor>a28905140f95bc587c3ffc03e253c108b</anchor>
      <arglist>(ocp_qp_condensing_config *config, void *dims_)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_condensing_calculate_size</name>
      <anchorfile>d6/dde/condensing__interface_8c.html</anchorfile>
      <anchor>a7d753bfb1cd9b092eeede9ea7ebc5fe3</anchor>
      <arglist>(ocp_qp_condensing_config *config, void *dims_, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>condensing_module *</type>
      <name>ocp_qp_condensing_assign</name>
      <anchorfile>d6/dde/condensing__interface_8c.html</anchorfile>
      <anchor>a1f1c13f1397f83655ee052c186c7ec37</anchor>
      <arglist>(ocp_qp_condensing_config *config, void *dims_, void *opts_, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>condensing_module *</type>
      <name>ocp_qp_condensing_create</name>
      <anchorfile>d6/dde/condensing__interface_8c.html</anchorfile>
      <anchor>a4b15b0fce448acb680eefcde3d6d2295</anchor>
      <arglist>(ocp_qp_condensing_config *config, void *dims_, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_condense</name>
      <anchorfile>d6/dde/condensing__interface_8c.html</anchorfile>
      <anchor>aaa4a67e1aef9540d5bb03f4e686fcf02</anchor>
      <arglist>(condensing_module *module, void *qp_in, void *qp_out)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_expand</name>
      <anchorfile>d6/dde/condensing__interface_8c.html</anchorfile>
      <anchor>a9826a871fd85f860a5ec89f7057ed309</anchor>
      <arglist>(condensing_module *module, void *qp_in, void *qp_out)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>condensing_interface.h</name>
    <path>/home/tom/sources/syscop/acados/interfaces/acados_c/</path>
    <filename>df/d4f/condensing__interface_8h</filename>
    <class kind="struct">condensing_plan</class>
    <class kind="struct">condensing_module</class>
    <member kind="enumeration">
      <type></type>
      <name>condensing_t</name>
      <anchorfile>df/d4f/condensing__interface_8h.html</anchorfile>
      <anchor>aaed6c542ba6fb47cf6a884084f8415c1</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>PARTIAL_CONDENSING</name>
      <anchorfile>df/d4f/condensing__interface_8h.html</anchorfile>
      <anchor>aaed6c542ba6fb47cf6a884084f8415c1a6c4013bc6af883c706e6de828ff2d6a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>FULL_CONDENSING</name>
      <anchorfile>df/d4f/condensing__interface_8h.html</anchorfile>
      <anchor>aaed6c542ba6fb47cf6a884084f8415c1a0837cb7ddde236a9d6867d55132dcbc1</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_condensing_config *</type>
      <name>ocp_qp_condensing_config_create</name>
      <anchorfile>df/d4f/condensing__interface_8h.html</anchorfile>
      <anchor>ae3a310a44fe882a61de503d27163c4c7</anchor>
      <arglist>(condensing_plan *plan)</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>ocp_qp_condensing_opts_create</name>
      <anchorfile>df/d4f/condensing__interface_8h.html</anchorfile>
      <anchor>a28905140f95bc587c3ffc03e253c108b</anchor>
      <arglist>(ocp_qp_condensing_config *config, void *dims_)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_condensing_calculate_size</name>
      <anchorfile>df/d4f/condensing__interface_8h.html</anchorfile>
      <anchor>a7d753bfb1cd9b092eeede9ea7ebc5fe3</anchor>
      <arglist>(ocp_qp_condensing_config *config, void *dims_, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>condensing_module *</type>
      <name>ocp_qp_condensing_assign</name>
      <anchorfile>df/d4f/condensing__interface_8h.html</anchorfile>
      <anchor>a1f1c13f1397f83655ee052c186c7ec37</anchor>
      <arglist>(ocp_qp_condensing_config *config, void *dims_, void *opts_, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>condensing_module *</type>
      <name>ocp_qp_condensing_create</name>
      <anchorfile>df/d4f/condensing__interface_8h.html</anchorfile>
      <anchor>a4b15b0fce448acb680eefcde3d6d2295</anchor>
      <arglist>(ocp_qp_condensing_config *config, void *dims_, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_condense</name>
      <anchorfile>df/d4f/condensing__interface_8h.html</anchorfile>
      <anchor>aaa4a67e1aef9540d5bb03f4e686fcf02</anchor>
      <arglist>(condensing_module *module, void *qp_in, void *qp_out)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_expand</name>
      <anchorfile>df/d4f/condensing__interface_8h.html</anchorfile>
      <anchor>a9826a871fd85f860a5ec89f7057ed309</anchor>
      <arglist>(condensing_module *module, void *qp_in, void *qp_out)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>dense_qp_interface.c</name>
    <path>/home/tom/sources/syscop/acados/interfaces/acados_c/</path>
    <filename>dc/d5e/dense__qp__interface_8c</filename>
    <includes id="dd/d7b/dense__qp__interface_8h" name="dense_qp_interface.h" local="yes" imported="no">acados_c/dense_qp_interface.h</includes>
    <member kind="function">
      <type>qp_solver_config *</type>
      <name>dense_qp_config_create</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>aa9feb39088298ba7ba10cd7abe20048a</anchor>
      <arglist>(dense_qp_solver_plan *plan)</arglist>
    </member>
    <member kind="function">
      <type>dense_qp_dims *</type>
      <name>dense_qp_dims_create</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>a02538dc6d2f50ddbceb7c0c6107db69e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>dense_qp_in *</type>
      <name>dense_qp_in_create</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>aa5a7795b607fc3aef8a7714c1fad1c85</anchor>
      <arglist>(qp_solver_config *config, dense_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>dense_qp_out *</type>
      <name>dense_qp_out_create</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>ad0eda6fd7d5096f064c2836f8d764d18</anchor>
      <arglist>(qp_solver_config *config, dense_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>dense_qp_opts_create</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>a04bad27017846dec1ccca3f1b309c41f</anchor>
      <arglist>(qp_solver_config *config, dense_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>dense_qp_calculate_size</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>a55b7018bf57cc170b378d44666657ee4</anchor>
      <arglist>(qp_solver_config *config, dense_qp_dims *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>dense_qp_solver *</type>
      <name>dense_qp_assign</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>a27057fc66c90ca8cd8b282176098e624</anchor>
      <arglist>(qp_solver_config *config, dense_qp_dims *dims, void *opts_, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>dense_qp_solver *</type>
      <name>dense_qp_create</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>ad2ee01f94ad051b49f76fc3a334bea4e</anchor>
      <arglist>(qp_solver_config *config, dense_qp_dims *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>dense_qp_solve</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>af9b8a60464c3fa91a43551132320b1d8</anchor>
      <arglist>(dense_qp_solver *solver, dense_qp_in *qp_in, dense_qp_out *qp_out)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static dense_qp_res *</type>
      <name>dense_qp_res_create</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>a2218fe8c782f9bca277e78a37c06708b</anchor>
      <arglist>(dense_qp_dims *dims)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static dense_qp_res_ws *</type>
      <name>dense_qp_res_workspace_create</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>a6e56af4ac7f111fa50adadfb6e65bdb0</anchor>
      <arglist>(dense_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>dense_qp_inf_norm_residuals</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>ac79abae6217985eb0cb1b547a86f8a36</anchor>
      <arglist>(dense_qp_dims *dims, dense_qp_in *qp_in, dense_qp_out *qp_out, double *res)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>dense_qp_set_field_double_array</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>a78a3a42655eab896c1213a4150a12d32</anchor>
      <arglist>(const char *field, double *arr, dense_qp_in *qp_in)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>dense_qp_set_field_int_array</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>aa7b78d10018097b8285297c1248e0998</anchor>
      <arglist>(const char *field, int *arr, dense_qp_in *qp_in)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>dense_qp_get_field_double_array</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>a183b7bdff570c05b8c14419066876655</anchor>
      <arglist>(const char *field, dense_qp_in *qp_in, double *arr)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>dense_qp_get_field_int_array</name>
      <anchorfile>dc/d5e/dense__qp__interface_8c.html</anchorfile>
      <anchor>a03f6547a4a68459fcb4273f6ce7e2c09</anchor>
      <arglist>(const char *field, dense_qp_in *qp_in, int *arr)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>dense_qp_interface.h</name>
    <path>/home/tom/sources/syscop/acados/interfaces/acados_c/</path>
    <filename>dd/d7b/dense__qp__interface_8h</filename>
    <includes id="d8/d1f/dense__qp__common_8h" name="dense_qp_common.h" local="yes" imported="no">acados/dense_qp/dense_qp_common.h</includes>
    <class kind="struct">dense_qp_solver_plan</class>
    <class kind="struct">dense_qp_solver</class>
    <member kind="enumeration">
      <type></type>
      <name>dense_qp_solver_t</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>a0c6d56ee417baf694b8709c0e7352604</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>DENSE_QP_HPIPM</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>a0c6d56ee417baf694b8709c0e7352604a8c8299b6d1ee3771855439d2f24b1294</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>DENSE_QP_QORE</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>a0c6d56ee417baf694b8709c0e7352604a4372a2dd17f3b71f358d12933f7f066e</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>DENSE_QP_QPOASES</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>a0c6d56ee417baf694b8709c0e7352604a7062c6652f6862f3916f0a11e7ff64d8</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>DENSE_QP_OOQP</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>a0c6d56ee417baf694b8709c0e7352604a3dd596cd00324886e85eca8d12550a88</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>qp_solver_config *</type>
      <name>dense_qp_config_create</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>aa9feb39088298ba7ba10cd7abe20048a</anchor>
      <arglist>(dense_qp_solver_plan *plan)</arglist>
    </member>
    <member kind="function">
      <type>dense_qp_dims *</type>
      <name>dense_qp_dims_create</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>a02538dc6d2f50ddbceb7c0c6107db69e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>dense_qp_in *</type>
      <name>dense_qp_in_create</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>aa5a7795b607fc3aef8a7714c1fad1c85</anchor>
      <arglist>(qp_solver_config *config, dense_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>dense_qp_out *</type>
      <name>dense_qp_out_create</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>ad0eda6fd7d5096f064c2836f8d764d18</anchor>
      <arglist>(qp_solver_config *config, dense_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>dense_qp_opts_create</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>a04bad27017846dec1ccca3f1b309c41f</anchor>
      <arglist>(qp_solver_config *config, dense_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>dense_qp_calculate_size</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>a55b7018bf57cc170b378d44666657ee4</anchor>
      <arglist>(qp_solver_config *config, dense_qp_dims *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>dense_qp_solver *</type>
      <name>dense_qp_assign</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>a27057fc66c90ca8cd8b282176098e624</anchor>
      <arglist>(qp_solver_config *config, dense_qp_dims *dims, void *opts_, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>dense_qp_solver *</type>
      <name>dense_qp_create</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>ad2ee01f94ad051b49f76fc3a334bea4e</anchor>
      <arglist>(qp_solver_config *config, dense_qp_dims *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>dense_qp_solve</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>af9b8a60464c3fa91a43551132320b1d8</anchor>
      <arglist>(dense_qp_solver *solver, dense_qp_in *qp_in, dense_qp_out *qp_out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>dense_qp_inf_norm_residuals</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>ac79abae6217985eb0cb1b547a86f8a36</anchor>
      <arglist>(dense_qp_dims *dims, dense_qp_in *qp_in, dense_qp_out *qp_out, double *res)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>dense_qp_set_field_double_array</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>a78a3a42655eab896c1213a4150a12d32</anchor>
      <arglist>(const char *field, double *arr, dense_qp_in *qp_in)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>dense_qp_set_field_int_array</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>aa7b78d10018097b8285297c1248e0998</anchor>
      <arglist>(const char *field, int *arr, dense_qp_in *qp_in)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>dense_qp_get_field_double_array</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>a183b7bdff570c05b8c14419066876655</anchor>
      <arglist>(const char *field, dense_qp_in *qp_in, double *arr)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>dense_qp_get_field_int_array</name>
      <anchorfile>dd/d7b/dense__qp__interface_8h.html</anchorfile>
      <anchor>a03f6547a4a68459fcb4273f6ce7e2c09</anchor>
      <arglist>(const char *field, dense_qp_in *qp_in, int *arr)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>external_function_interface.c</name>
    <path>/home/tom/sources/syscop/acados/interfaces/acados_c/</path>
    <filename>d3/d8c/external__function__interface_8c</filename>
    <includes id="d2/d7a/external__function__interface_8h" name="external_function_interface.h" local="yes" imported="no">acados_c/external_function_interface.h</includes>
    <includes id="dd/ddf/external__function__generic_8h" name="external_function_generic.h" local="yes" imported="no">acados/utils/external_function_generic.h</includes>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_create</name>
      <anchorfile>d3/d8c/external__function__interface_8c.html</anchorfile>
      <anchor>a71db68392aa9fd6ac333b97bc8abc9e9</anchor>
      <arglist>(external_function_casadi *fun)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_create_array</name>
      <anchorfile>d3/d8c/external__function__interface_8c.html</anchorfile>
      <anchor>a1472bd8b31cdd1b09ea7157b8707d862</anchor>
      <arglist>(int size, external_function_casadi *funs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_free</name>
      <anchorfile>d3/d8c/external__function__interface_8c.html</anchorfile>
      <anchor>a1b1dafbb6cc975c8d9efda8b7d3f0773</anchor>
      <arglist>(external_function_casadi *fun)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_free_array</name>
      <anchorfile>d3/d8c/external__function__interface_8c.html</anchorfile>
      <anchor>af1219c0773ccb3c5d0fc7527d7034b3a</anchor>
      <arglist>(int size, external_function_casadi *funs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_create</name>
      <anchorfile>d3/d8c/external__function__interface_8c.html</anchorfile>
      <anchor>acf71abfb32102b70a5bea140dddb7086</anchor>
      <arglist>(external_function_param_casadi *fun, int np)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_create_array</name>
      <anchorfile>d3/d8c/external__function__interface_8c.html</anchorfile>
      <anchor>ada7fece668ebec795383a4c2ff25471d</anchor>
      <arglist>(int size, external_function_param_casadi *funs, int np)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_free</name>
      <anchorfile>d3/d8c/external__function__interface_8c.html</anchorfile>
      <anchor>a8758793f59675efd4881824bf727d421</anchor>
      <arglist>(external_function_param_casadi *fun)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_free_array</name>
      <anchorfile>d3/d8c/external__function__interface_8c.html</anchorfile>
      <anchor>a2eb8cf6454ba161d11fff4f001b5d005</anchor>
      <arglist>(int size, external_function_param_casadi *funs)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>external_function_interface.h</name>
    <path>/home/tom/sources/syscop/acados/interfaces/acados_c/</path>
    <filename>d2/d7a/external__function__interface_8h</filename>
    <includes id="dd/ddf/external__function__generic_8h" name="external_function_generic.h" local="yes" imported="no">acados/utils/external_function_generic.h</includes>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_create</name>
      <anchorfile>d2/d7a/external__function__interface_8h.html</anchorfile>
      <anchor>a71db68392aa9fd6ac333b97bc8abc9e9</anchor>
      <arglist>(external_function_casadi *fun)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_free</name>
      <anchorfile>d2/d7a/external__function__interface_8h.html</anchorfile>
      <anchor>a1b1dafbb6cc975c8d9efda8b7d3f0773</anchor>
      <arglist>(external_function_casadi *fun)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_create_array</name>
      <anchorfile>d2/d7a/external__function__interface_8h.html</anchorfile>
      <anchor>a1472bd8b31cdd1b09ea7157b8707d862</anchor>
      <arglist>(int size, external_function_casadi *funs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_casadi_free_array</name>
      <anchorfile>d2/d7a/external__function__interface_8h.html</anchorfile>
      <anchor>af1219c0773ccb3c5d0fc7527d7034b3a</anchor>
      <arglist>(int size, external_function_casadi *funs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_create</name>
      <anchorfile>d2/d7a/external__function__interface_8h.html</anchorfile>
      <anchor>acf71abfb32102b70a5bea140dddb7086</anchor>
      <arglist>(external_function_param_casadi *fun, int np)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_free</name>
      <anchorfile>d2/d7a/external__function__interface_8h.html</anchorfile>
      <anchor>a8758793f59675efd4881824bf727d421</anchor>
      <arglist>(external_function_param_casadi *fun)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_create_array</name>
      <anchorfile>d2/d7a/external__function__interface_8h.html</anchorfile>
      <anchor>ada7fece668ebec795383a4c2ff25471d</anchor>
      <arglist>(int size, external_function_param_casadi *funs, int np)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>external_function_param_casadi_free_array</name>
      <anchorfile>d2/d7a/external__function__interface_8h.html</anchorfile>
      <anchor>a2eb8cf6454ba161d11fff4f001b5d005</anchor>
      <arglist>(int size, external_function_param_casadi *funs)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>ocp_nlp_interface.c</name>
    <path>/home/tom/sources/syscop/acados/interfaces/acados_c/</path>
    <filename>db/d2b/ocp__nlp__interface_8c</filename>
    <includes id="d3/df7/ocp__nlp__interface_8h" name="ocp_nlp_interface.h" local="yes" imported="no">acados_c/ocp_nlp_interface.h</includes>
    <member kind="function" static="yes">
      <type>static int</type>
      <name>ocp_nlp_plan_calculate_size</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a53e099ad51079b0a23b47d8ee1065661</anchor>
      <arglist>(int N)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static ocp_nlp_plan *</type>
      <name>ocp_nlp_plan_assign</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>ab3e53ff141952734431f83d93e56ca3f</anchor>
      <arglist>(int N, void *raw_memory)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static void</type>
      <name>ocp_nlp_plan_initialize_default</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a6cb64eda98cd1f272ff25bbe4f7f3b3c</anchor>
      <arglist>(ocp_nlp_plan *plan)</arglist>
    </member>
    <member kind="function">
      <type>ocp_nlp_plan *</type>
      <name>ocp_nlp_plan_create</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a4a2848bc3832f1813c1cc36a6a83f3da</anchor>
      <arglist>(int N)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_plan_destroy</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a67269f69bd291c8a3da63f4fc46d8ed9</anchor>
      <arglist>(void *plan_)</arglist>
    </member>
    <member kind="function">
      <type>ocp_nlp_config *</type>
      <name>ocp_nlp_config_create</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>af079733fd52ec57419778eacda6deda6</anchor>
      <arglist>(ocp_nlp_plan plan)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_config_destroy</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a89ebfa6845c7c4bf431fd367d13fc300</anchor>
      <arglist>(void *config_)</arglist>
    </member>
    <member kind="function">
      <type>ocp_nlp_dims *</type>
      <name>ocp_nlp_dims_create</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>ab8ba2218add35550ec3afdbfc359a096</anchor>
      <arglist>(void *config_)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_dims_destroy</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a66baa50d39d59f8dcaa3f2bbc2658ca6</anchor>
      <arglist>(void *dims_)</arglist>
    </member>
    <member kind="function">
      <type>ocp_nlp_in *</type>
      <name>ocp_nlp_in_create</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a76c45dc6a07eeaee7b012fa2720b2961</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_in_destroy</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>aa7e9da740e16a057e74646e9dd72bd03</anchor>
      <arglist>(void *in)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_in_set</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a5b9b022caedf8903bebafebc98e13c77</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, int stage, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_nlp_dynamics_model_set</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a8fce965b258aeeea550d07e028320db2</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, int stage, const char *fun_type, void *fun_ptr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_nlp_cost_model_set</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>abdff6e8c13af265aaab3c62dbb2667d8</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, int stage, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>nlp_set_discrete_model_in_stage</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>ab1b555c8260df8766ed7c783a3c1650d</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_in *in, int stage, void *fun_ptr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_nlp_constraints_model_set</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a1e57bc3b81907c9015ad8a60838aeac6</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, int stage, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>ocp_nlp_out *</type>
      <name>ocp_nlp_out_create</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a0acab204cf41287c8cce5230a57a5eca</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_out_destroy</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>ab964e9550c2d6d9e4d99e650d1b2ca53</anchor>
      <arglist>(void *out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_out_set</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a611fc758dc86424dde430f8a49d1c372</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out, int stage, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_out_get</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a2d3b7666a1907c6f8454580274a332ff</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out, int stage, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>ocp_nlp_opts_create</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a1cffcd139bb98d7d2f56004c57f588e9</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_opts_set</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>ab185e82989f3f09378863050fb17561a</anchor>
      <arglist>(ocp_nlp_config *config, void *opts_, const char *field, const void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_dynamics_opts_set</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>aecdd911a7739308f097266d93208e19c</anchor>
      <arglist>(ocp_nlp_config *config, void *opts_, int stage, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_opts_update</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a2152c81d02c80272c5113445181b7eab</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_opts_destroy</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>aea5529907e3348acc38506ead8e2f6f3</anchor>
      <arglist>(void *opts)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static int</type>
      <name>ocp_nlp_calculate_size</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a98d5365e36a52fef6f4fa5e71a5d857a</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static ocp_nlp_solver *</type>
      <name>ocp_nlp_assign</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a592c4c481f50a04ac0d61cb375327d0b</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>ocp_nlp_solver *</type>
      <name>ocp_nlp_solver_create</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a070da74bbc8ddf921b82ae22719084c4</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_solver_destroy</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>ab83da7685e510c43e0efe6769d349eec</anchor>
      <arglist>(void *solver)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_nlp_solve</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a1fac9ed35c5282447fc80c12d74519a7</anchor>
      <arglist>(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_nlp_precompute</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>aa3ac6336a05ef30bf5f0f1c536273727</anchor>
      <arglist>(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_get</name>
      <anchorfile>db/d2b/ocp__nlp__interface_8c.html</anchorfile>
      <anchor>a9eea79aa2137d46bf33a490172d00d4f</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_solver *solver, const char *field, void *return_value_)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>ocp_nlp_interface.h</name>
    <path>/home/tom/sources/syscop/acados/interfaces/acados_c/</path>
    <filename>d3/df7/ocp__nlp__interface_8h</filename>
    <includes id="d4/db9/ocp__qp__interface_8h" name="ocp_qp_interface.h" local="yes" imported="no">acados_c/ocp_qp_interface.h</includes>
    <includes id="d3/d69/sim__interface_8h" name="sim_interface.h" local="yes" imported="no">acados_c/sim_interface.h</includes>
    <class kind="struct">ocp_nlp_plan</class>
    <class kind="struct">ocp_nlp_solver</class>
    <member kind="enumeration">
      <type></type>
      <name>ocp_nlp_solver_t</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>ab8163075e7393832509b38f3f5a3095b</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>SQP</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>ab8163075e7393832509b38f3f5a3095baca52627dda20bfd9180bdad14194e39a</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>SQP_RTI</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>ab8163075e7393832509b38f3f5a3095bafd9c575a106d3f88c997a2fcc7677b96</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>INVALID_NLP_SOLVER</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>ab8163075e7393832509b38f3f5a3095badae9e964c4c5073488a5b231a3c014d4</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumeration">
      <type></type>
      <name>ocp_nlp_cost_t</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a7758d4da081915969df6996e1208b52e</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>LINEAR_LS</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a7758d4da081915969df6996e1208b52ea6d5de376d51f43b79dcde02a807d6432</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>NONLINEAR_LS</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a7758d4da081915969df6996e1208b52ea977ff058ad9f2008b45f96b0af38b225</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>EXTERNALLY_PROVIDED</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a7758d4da081915969df6996e1208b52ea8f3b9cf836b0048b18c0b71020df5254</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>INVALID_COST</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a7758d4da081915969df6996e1208b52ea966cca8bc1628828a945af66d94f77b5</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumeration">
      <type></type>
      <name>ocp_nlp_dynamics_t</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a77b6b533423c73460e54008f9d342f3a</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>CONTINUOUS_MODEL</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a77b6b533423c73460e54008f9d342f3aac29ab7ad29ea7fed429e8ab76d99951a</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>DISCRETE_MODEL</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a77b6b533423c73460e54008f9d342f3aaf630e1f3bf526545d6d8378abc5dddc1</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>INVALID_DYNAMICS</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a77b6b533423c73460e54008f9d342f3aa2ffac5d3bdf3282d38251b5104579172</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumeration">
      <type></type>
      <name>ocp_nlp_constraints_t</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a798ddd04e62b1b6534dd957c35a0cb9e</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>BGH</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a798ddd04e62b1b6534dd957c35a0cb9ea1d3b30250442be5cf6d0a6f243ebd6d0</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>BGHP</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a798ddd04e62b1b6534dd957c35a0cb9eae9fe10129638e3a3f28e587ee823fab6</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>INVALID_CONSTRAINT</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a798ddd04e62b1b6534dd957c35a0cb9ea25c1070d6690ff92e17b6e2c69c2e34c</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumeration">
      <type></type>
      <name>ocp_nlp_reg_t</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>ad6df73d48744b682a7853f9a6a21335c</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>NO_REGULARIZE</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>ad6df73d48744b682a7853f9a6a21335ca5b5cc5b4fe50f3dc084aa4e2c476328b</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>MIRROR</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>ad6df73d48744b682a7853f9a6a21335cab9aa0079c7544a4246b445bb04ef0352</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>PROJECT</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>ad6df73d48744b682a7853f9a6a21335ca27d5e6e2e38089805b568b424d86587d</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>CONVEXIFY</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>ad6df73d48744b682a7853f9a6a21335ca2230e4bc3f6584d207c12d644299f8dc</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>INVALID_REGULARIZE</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>ad6df73d48744b682a7853f9a6a21335cab3f4ab052e17292f80ced157fdde9bbe</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>ocp_nlp_plan *</type>
      <name>ocp_nlp_plan_create</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a4a2848bc3832f1813c1cc36a6a83f3da</anchor>
      <arglist>(int N)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_plan_destroy</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a67269f69bd291c8a3da63f4fc46d8ed9</anchor>
      <arglist>(void *plan_)</arglist>
    </member>
    <member kind="function">
      <type>ocp_nlp_config *</type>
      <name>ocp_nlp_config_create</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>af079733fd52ec57419778eacda6deda6</anchor>
      <arglist>(ocp_nlp_plan plan)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_config_destroy</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a89ebfa6845c7c4bf431fd367d13fc300</anchor>
      <arglist>(void *config_)</arglist>
    </member>
    <member kind="function">
      <type>ocp_nlp_dims *</type>
      <name>ocp_nlp_dims_create</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>ab8ba2218add35550ec3afdbfc359a096</anchor>
      <arglist>(void *config_)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_dims_destroy</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a66baa50d39d59f8dcaa3f2bbc2658ca6</anchor>
      <arglist>(void *dims_)</arglist>
    </member>
    <member kind="function">
      <type>ocp_nlp_in *</type>
      <name>ocp_nlp_in_create</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a76c45dc6a07eeaee7b012fa2720b2961</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_in_destroy</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>aa7e9da740e16a057e74646e9dd72bd03</anchor>
      <arglist>(void *in)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_in_set</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a5b9b022caedf8903bebafebc98e13c77</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, int stage, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_nlp_dynamics_model_set</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a8fce965b258aeeea550d07e028320db2</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, int stage, const char *fun_type, void *fun_ptr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>nlp_set_discrete_model_in_stage</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>ab1b555c8260df8766ed7c783a3c1650d</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_in *in, int stage, void *fun_ptr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_nlp_cost_model_set</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>abdff6e8c13af265aaab3c62dbb2667d8</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, int stage, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_nlp_constraints_model_set</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a1e57bc3b81907c9015ad8a60838aeac6</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, int stage, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>ocp_nlp_out *</type>
      <name>ocp_nlp_out_create</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a0acab204cf41287c8cce5230a57a5eca</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_out_destroy</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>ab964e9550c2d6d9e4d99e650d1b2ca53</anchor>
      <arglist>(void *out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_out_set</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a611fc758dc86424dde430f8a49d1c372</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out, int stage, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_out_get</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a2d3b7666a1907c6f8454580274a332ff</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out, int stage, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>ocp_nlp_opts_create</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a1cffcd139bb98d7d2f56004c57f588e9</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_opts_destroy</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>aea5529907e3348acc38506ead8e2f6f3</anchor>
      <arglist>(void *opts)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_opts_set</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>ab185e82989f3f09378863050fb17561a</anchor>
      <arglist>(ocp_nlp_config *config, void *opts_, const char *field, const void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_dynamics_opts_set</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>aecdd911a7739308f097266d93208e19c</anchor>
      <arglist>(ocp_nlp_config *config, void *opts_, int stage, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_opts_update</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a2152c81d02c80272c5113445181b7eab</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>ocp_nlp_solver *</type>
      <name>ocp_nlp_solver_create</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a070da74bbc8ddf921b82ae22719084c4</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_solver_destroy</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>ab83da7685e510c43e0efe6769d349eec</anchor>
      <arglist>(void *solver)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_nlp_solve</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a1fac9ed35c5282447fc80c12d74519a7</anchor>
      <arglist>(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_nlp_precompute</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>aa3ac6336a05ef30bf5f0f1c536273727</anchor>
      <arglist>(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_nlp_get</name>
      <anchorfile>d3/df7/ocp__nlp__interface_8h.html</anchorfile>
      <anchor>a9eea79aa2137d46bf33a490172d00d4f</anchor>
      <arglist>(ocp_nlp_config *config, ocp_nlp_solver *solver, const char *field, void *return_value_)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>ocp_qp_interface.c</name>
    <path>/home/tom/sources/syscop/acados/interfaces/acados_c/</path>
    <filename>d3/d86/ocp__qp__interface_8c</filename>
    <includes id="d4/db9/ocp__qp__interface_8h" name="ocp_qp_interface.h" local="yes" imported="no">acados_c/ocp_qp_interface.h</includes>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_xcond_solver_config_initialize_default</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>a4a9a32dde33de19f904c6105b7757aaa</anchor>
      <arglist>(ocp_qp_solver_t solver_name, ocp_qp_xcond_solver_config *solver_config)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_xcond_solver_config *</type>
      <name>ocp_qp_config_create</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>ad0c4fee21f57f33429283352c0fdc23b</anchor>
      <arglist>(ocp_qp_solver_plan plan)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_config_free</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>a57831bb14422c431c9cfd19ad60e98c7</anchor>
      <arglist>(void *config)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_dims *</type>
      <name>ocp_qp_dims_create</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>a351dc1d43c859234277bb269bdb8706f</anchor>
      <arglist>(int N)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_dims_free</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>ab8f90af8d06f2a38a236b2554cbfde72</anchor>
      <arglist>(void *dims_)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_in *</type>
      <name>ocp_qp_in_create</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>a447809cc06aba6f4be49973ef7b34877</anchor>
      <arglist>(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_in_free</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>a583a1f198497cb3134970f0278e580ff</anchor>
      <arglist>(void *in_)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_out *</type>
      <name>ocp_qp_out_create</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>a8ffd49cadb2aad910cc5a8a166a4043c</anchor>
      <arglist>(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_out_free</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>a0f49bd988087d3c6b2996b7f2a8bd6cc</anchor>
      <arglist>(void *out_)</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>ocp_qp_opts_create</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>a7470b0e6eef9dde0cab76576ea715ee4</anchor>
      <arglist>(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_opts_free</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>a21880271a6758c276e9d614d8a11ebb5</anchor>
      <arglist>(void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_calculate_size</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>ae6db957993cc71c3ca1f79a97d84dcb1</anchor>
      <arglist>(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_solver *</type>
      <name>ocp_qp_assign</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>a50d96db297aba82acded233ea5d8dc8e</anchor>
      <arglist>(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims, void *opts_, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_solver *</type>
      <name>ocp_qp_create</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>ad3f84b7a0ea53e54657da9407aa1d636</anchor>
      <arglist>(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_solve</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>af65c078579ef56546a5284ebb537ee29</anchor>
      <arglist>(ocp_qp_solver *solver, ocp_qp_in *qp_in, ocp_qp_out *qp_out)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static ocp_qp_res *</type>
      <name>ocp_qp_res_create</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>aaf69f8c623d0cd85d31eb1fe7e9f3308</anchor>
      <arglist>(ocp_qp_dims *dims)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static ocp_qp_res_ws *</type>
      <name>ocp_qp_res_workspace_create</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>ad13cd631147a4927b192221a39a01e31</anchor>
      <arglist>(ocp_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_inf_norm_residuals</name>
      <anchorfile>d3/d86/ocp__qp__interface_8c.html</anchorfile>
      <anchor>aa4a8b143fda7f600bbeebce13d063487</anchor>
      <arglist>(ocp_qp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out, double *res)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>ocp_qp_interface.h</name>
    <path>/home/tom/sources/syscop/acados/interfaces/acados_c/</path>
    <filename>d4/db9/ocp__qp__interface_8h</filename>
    <includes id="d5/d0d/ocp__qp__common_8h" name="ocp_qp_common.h" local="yes" imported="no">acados/ocp_qp/ocp_qp_common.h</includes>
    <class kind="struct">ocp_qp_solver_plan</class>
    <class kind="struct">ocp_qp_solver</class>
    <member kind="enumeration">
      <type></type>
      <name>ocp_qp_solver_t</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>a9ef89456f459c1d044afbf03b560e950</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>PARTIAL_CONDENSING_HPIPM</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>a9ef89456f459c1d044afbf03b560e950a7baa9d9ce9bdf39b38a831e327f7ee67</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>FULL_CONDENSING_HPIPM</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>a9ef89456f459c1d044afbf03b560e950a36d6664f93a24f14cb86b9d64b85d378</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>INVALID_QP_SOLVER</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>a9ef89456f459c1d044afbf03b560e950a32c24a4f759f6c18065cbabd0c75bd6d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_xcond_solver_config_initialize_default</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>a4a9a32dde33de19f904c6105b7757aaa</anchor>
      <arglist>(ocp_qp_solver_t solver_name, ocp_qp_xcond_solver_config *solver_config)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_xcond_solver_config *</type>
      <name>ocp_qp_config_create</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>ad0c4fee21f57f33429283352c0fdc23b</anchor>
      <arglist>(ocp_qp_solver_plan plan)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_config_free</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>a8f9d66f7b3cbe1db5444ac0befe27b5f</anchor>
      <arglist>(void *config_)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_dims *</type>
      <name>ocp_qp_dims_create</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>a351dc1d43c859234277bb269bdb8706f</anchor>
      <arglist>(int N)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_dims_free</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>ab8f90af8d06f2a38a236b2554cbfde72</anchor>
      <arglist>(void *dims_)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_in *</type>
      <name>ocp_qp_in_create</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>a447809cc06aba6f4be49973ef7b34877</anchor>
      <arglist>(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_in_free</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>a583a1f198497cb3134970f0278e580ff</anchor>
      <arglist>(void *in_)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_out *</type>
      <name>ocp_qp_out_create</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>a8ffd49cadb2aad910cc5a8a166a4043c</anchor>
      <arglist>(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_out_free</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>a0f49bd988087d3c6b2996b7f2a8bd6cc</anchor>
      <arglist>(void *out_)</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>ocp_qp_opts_create</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>a7470b0e6eef9dde0cab76576ea715ee4</anchor>
      <arglist>(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_opts_free</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>a21880271a6758c276e9d614d8a11ebb5</anchor>
      <arglist>(void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_calculate_size</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>ae6db957993cc71c3ca1f79a97d84dcb1</anchor>
      <arglist>(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_solver *</type>
      <name>ocp_qp_assign</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>a50d96db297aba82acded233ea5d8dc8e</anchor>
      <arglist>(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims, void *opts_, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>ocp_qp_solver *</type>
      <name>ocp_qp_create</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>ad3f84b7a0ea53e54657da9407aa1d636</anchor>
      <arglist>(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ocp_qp_solve</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>af65c078579ef56546a5284ebb537ee29</anchor>
      <arglist>(ocp_qp_solver *solver, ocp_qp_in *qp_in, ocp_qp_out *qp_out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ocp_qp_inf_norm_residuals</name>
      <anchorfile>d4/db9/ocp__qp__interface_8h.html</anchorfile>
      <anchor>aa4a8b143fda7f600bbeebce13d063487</anchor>
      <arglist>(ocp_qp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out, double *res)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>options_interface.c</name>
    <path>/home/tom/sources/syscop/acados/interfaces/acados_c/</path>
    <filename>d7/d4f/options__interface_8c</filename>
    <includes id="dc/dac/options__interface_8h" name="options_interface.h" local="yes" imported="no">acados_c/options_interface.h</includes>
    <member kind="function" static="yes">
      <type>static bool</type>
      <name>is_qp_solver</name>
      <anchorfile>d7/d4f/options__interface_8c.html</anchorfile>
      <anchor>aff2baf9b79205d728458bec8497e610e</anchor>
      <arglist>(const char *qp_solver_name)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_option_int</name>
      <anchorfile>d7/d4f/options__interface_8c.html</anchorfile>
      <anchor>ae252413931289ccd9803cbf5536b2404</anchor>
      <arglist>(const void *args_, const char *option)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>set_option_int</name>
      <anchorfile>d7/d4f/options__interface_8c.html</anchorfile>
      <anchor>a65b808d8b30ea5aeb835efbf51381c0f</anchor>
      <arglist>(void *args_, const char *option, const int value)</arglist>
    </member>
    <member kind="function">
      <type>const int *</type>
      <name>get_option_int_array</name>
      <anchorfile>d7/d4f/options__interface_8c.html</anchorfile>
      <anchor>ae0cfed1ac1148e2a2dd8f88728ff44df</anchor>
      <arglist>(const void *args_, const char *option)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>set_option_int_array</name>
      <anchorfile>d7/d4f/options__interface_8c.html</anchorfile>
      <anchor>a353febe263cd08b748ec7a07b176f1de</anchor>
      <arglist>(void *args_, const char *option, const int *value)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_option_double</name>
      <anchorfile>d7/d4f/options__interface_8c.html</anchorfile>
      <anchor>a7bdb6122b90762da3e98b3be9dd9a2f9</anchor>
      <arglist>(const void *args_, const char *option)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>set_option_double</name>
      <anchorfile>d7/d4f/options__interface_8c.html</anchorfile>
      <anchor>a0443644a9a035b927da640e58925eec8</anchor>
      <arglist>(void *args_, const char *option, const double value)</arglist>
    </member>
    <member kind="function">
      <type>const double *</type>
      <name>get_option_double_array</name>
      <anchorfile>d7/d4f/options__interface_8c.html</anchorfile>
      <anchor>aa26fe63c9b7188c5286296e1fe8bc6b3</anchor>
      <arglist>(const void *args_, const char *option)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>set_option_double_array</name>
      <anchorfile>d7/d4f/options__interface_8c.html</anchorfile>
      <anchor>af6a1310919dd7f40c435f3badf797932</anchor>
      <arglist>(void *args_, const char *option, const double *value)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>options_interface.h</name>
    <path>/home/tom/sources/syscop/acados/interfaces/acados_c/</path>
    <filename>dc/dac/options__interface_8h</filename>
    <member kind="function">
      <type>int</type>
      <name>get_option_int</name>
      <anchorfile>dc/dac/options__interface_8h.html</anchorfile>
      <anchor>ae252413931289ccd9803cbf5536b2404</anchor>
      <arglist>(const void *args_, const char *option)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>set_option_int</name>
      <anchorfile>dc/dac/options__interface_8h.html</anchorfile>
      <anchor>a65b808d8b30ea5aeb835efbf51381c0f</anchor>
      <arglist>(void *args_, const char *option, const int value)</arglist>
    </member>
    <member kind="function">
      <type>const int *</type>
      <name>get_option_int_array</name>
      <anchorfile>dc/dac/options__interface_8h.html</anchorfile>
      <anchor>ae0cfed1ac1148e2a2dd8f88728ff44df</anchor>
      <arglist>(const void *args_, const char *option)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>set_option_int_array</name>
      <anchorfile>dc/dac/options__interface_8h.html</anchorfile>
      <anchor>a353febe263cd08b748ec7a07b176f1de</anchor>
      <arglist>(void *args_, const char *option, const int *value)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_option_double</name>
      <anchorfile>dc/dac/options__interface_8h.html</anchorfile>
      <anchor>a7bdb6122b90762da3e98b3be9dd9a2f9</anchor>
      <arglist>(const void *args_, const char *option)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>set_option_double</name>
      <anchorfile>dc/dac/options__interface_8h.html</anchorfile>
      <anchor>a0443644a9a035b927da640e58925eec8</anchor>
      <arglist>(void *args_, const char *option, const double value)</arglist>
    </member>
    <member kind="function">
      <type>const double *</type>
      <name>get_option_double_array</name>
      <anchorfile>dc/dac/options__interface_8h.html</anchorfile>
      <anchor>aa26fe63c9b7188c5286296e1fe8bc6b3</anchor>
      <arglist>(const void *args_, const char *option)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>set_option_double_array</name>
      <anchorfile>dc/dac/options__interface_8h.html</anchorfile>
      <anchor>af6a1310919dd7f40c435f3badf797932</anchor>
      <arglist>(void *args_, const char *option, const double *value)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>sim_interface.c</name>
    <path>/home/tom/sources/syscop/acados/interfaces/acados_c/</path>
    <filename>d8/d97/sim__interface_8c</filename>
    <includes id="dc/d3b/sim__common_8h" name="sim_common.h" local="yes" imported="no">acados/sim/sim_common.h</includes>
    <includes id="d3/d69/sim__interface_8h" name="sim_interface.h" local="yes" imported="no">acados_c/sim_interface.h</includes>
    <member kind="function">
      <type>sim_config *</type>
      <name>sim_config_create</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>ac6815b30545301745714353a8cfce71d</anchor>
      <arglist>(sim_solver_plan plan)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_config_destroy</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>a284ce61a2ba160f9cd3e5be3c44d416c</anchor>
      <arglist>(void *config)</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>sim_dims_create</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>a56e49efaf75fc3e429f67477bdf9bacc</anchor>
      <arglist>(void *config_)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_dims_destroy</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>ae113d216119b8170b9c2d1323256a5cc</anchor>
      <arglist>(void *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_dims_set</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>a73296134b77fb9735df1cf8d5f91fe80</anchor>
      <arglist>(sim_config *config, void *dims, const char *field, const int *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_dims_get</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>a3001c25f6e1e5fb6f73705f434f1ec37</anchor>
      <arglist>(sim_config *config, void *dims, const char *field, int *value)</arglist>
    </member>
    <member kind="function">
      <type>sim_in *</type>
      <name>sim_in_create</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>a36b505391348e85915faba8111553922</anchor>
      <arglist>(sim_config *config, void *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_in_destroy</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>a7388419838155d6ff7ae16b9d4312284</anchor>
      <arglist>(void *in)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_in_set</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>ad58ac86a6346973e313eb515cc1c12e7</anchor>
      <arglist>(void *config_, void *dims_, sim_in *in, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>sim_out *</type>
      <name>sim_out_create</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>aae2cf940a861afb282121781a4c15376</anchor>
      <arglist>(sim_config *config, void *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_out_destroy</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>ad7091abb99af135a01c753485129e89f</anchor>
      <arglist>(void *out)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_out_get</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>a455052a9d56137d4a2ae7668126e6067</anchor>
      <arglist>(void *config, void *dims, sim_out *out, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>sim_opts_create</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>a19ab776ed52ab128a48e22557397f368</anchor>
      <arglist>(sim_config *config, void *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_opts_destroy</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>abeb24edd4a5c03ec121cb0dfaac1a1e2</anchor>
      <arglist>(void *opts)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_opts_set</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>a32dd1e905140a0a79c953474ec2a22e2</anchor>
      <arglist>(sim_config *config, void *opts, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_calculate_size</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>a669c972999a68e1d397a657d13d8de89</anchor>
      <arglist>(sim_config *config, void *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>sim_solver *</type>
      <name>sim_assign</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>a5d8d673e138d0543fbb7b98047d1f9dc</anchor>
      <arglist>(sim_config *config, void *dims, void *opts_, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>sim_solver *</type>
      <name>sim_solver_create</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>a4bfdf874c97ed6a9349ba4a4412c3cb6</anchor>
      <arglist>(sim_config *config, void *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_solver_destroy</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>af1e4270ed72f0023438924ffa75ce1f4</anchor>
      <arglist>(void *solver)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_solve</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>a16b553e998f20c93119748454c533505</anchor>
      <arglist>(sim_solver *solver, sim_in *in, sim_out *out)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_precompute</name>
      <anchorfile>d8/d97/sim__interface_8c.html</anchorfile>
      <anchor>aa8f55d61e23edd951db0435f7c49f692</anchor>
      <arglist>(sim_solver *solver, sim_in *in, sim_out *out)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>sim_interface.h</name>
    <path>/home/tom/sources/syscop/acados/interfaces/acados_c/</path>
    <filename>d3/d69/sim__interface_8h</filename>
    <includes id="dc/d3b/sim__common_8h" name="sim_common.h" local="yes" imported="no">acados/sim/sim_common.h</includes>
    <class kind="struct">sim_solver_plan</class>
    <class kind="struct">sim_solver</class>
    <member kind="enumeration">
      <type></type>
      <name>sim_solver_t</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>ac280abde45c53996ed0ef65412e2ff1a</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>ERK</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>ac280abde45c53996ed0ef65412e2ff1aa16e08b722a8d1f3ca12a88098b959326</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>IRK</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>ac280abde45c53996ed0ef65412e2ff1aabefbc9c63d7cd5150fa7d7a4fabd30cd</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>GNSF</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>ac280abde45c53996ed0ef65412e2ff1aa9d78424a22d43f540bf2958e2b5cf04d</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>LIFTED_IRK</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>ac280abde45c53996ed0ef65412e2ff1aa3756bfcabeb364e52ad47bc7fe760bf3</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>INVALID_SIM_SOLVER</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>ac280abde45c53996ed0ef65412e2ff1aaadecbaceb6f49e7636be6bb1b5855efb</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>sim_config *</type>
      <name>sim_config_create</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>ac6815b30545301745714353a8cfce71d</anchor>
      <arglist>(sim_solver_plan plan)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_config_destroy</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>a284ce61a2ba160f9cd3e5be3c44d416c</anchor>
      <arglist>(void *config)</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>sim_dims_create</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>a56e49efaf75fc3e429f67477bdf9bacc</anchor>
      <arglist>(void *config_)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_dims_destroy</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>ae113d216119b8170b9c2d1323256a5cc</anchor>
      <arglist>(void *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_dims_set</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>a73296134b77fb9735df1cf8d5f91fe80</anchor>
      <arglist>(sim_config *config, void *dims, const char *field, const int *value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_dims_get</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>a3001c25f6e1e5fb6f73705f434f1ec37</anchor>
      <arglist>(sim_config *config, void *dims, const char *field, int *value)</arglist>
    </member>
    <member kind="function">
      <type>sim_in *</type>
      <name>sim_in_create</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>a36b505391348e85915faba8111553922</anchor>
      <arglist>(sim_config *config, void *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_in_destroy</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>aa1007ce29e9b89376db4f5292bba09d3</anchor>
      <arglist>(void *out)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_in_set</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>ad58ac86a6346973e313eb515cc1c12e7</anchor>
      <arglist>(void *config_, void *dims_, sim_in *in, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>sim_out *</type>
      <name>sim_out_create</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>aae2cf940a861afb282121781a4c15376</anchor>
      <arglist>(sim_config *config, void *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_out_destroy</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>ad7091abb99af135a01c753485129e89f</anchor>
      <arglist>(void *out)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_out_get</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>a455052a9d56137d4a2ae7668126e6067</anchor>
      <arglist>(void *config, void *dims, sim_out *out, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>sim_opts_create</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>a19ab776ed52ab128a48e22557397f368</anchor>
      <arglist>(sim_config *config, void *dims)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_opts_destroy</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>abeb24edd4a5c03ec121cb0dfaac1a1e2</anchor>
      <arglist>(void *opts)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_opts_set</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>a32dd1e905140a0a79c953474ec2a22e2</anchor>
      <arglist>(sim_config *config, void *opts, const char *field, void *value)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_calculate_size</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>a669c972999a68e1d397a657d13d8de89</anchor>
      <arglist>(sim_config *config, void *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>sim_solver *</type>
      <name>sim_assign</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>a5d8d673e138d0543fbb7b98047d1f9dc</anchor>
      <arglist>(sim_config *config, void *dims, void *opts_, void *raw_memory)</arglist>
    </member>
    <member kind="function">
      <type>sim_solver *</type>
      <name>sim_solver_create</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>a4bfdf874c97ed6a9349ba4a4412c3cb6</anchor>
      <arglist>(sim_config *config, void *dims, void *opts_)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>sim_solver_destroy</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>af1e4270ed72f0023438924ffa75ce1f4</anchor>
      <arglist>(void *solver)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_solve</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>a16b553e998f20c93119748454c533505</anchor>
      <arglist>(sim_solver *solver, sim_in *in, sim_out *out)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sim_precompute</name>
      <anchorfile>d3/d69/sim__interface_8h.html</anchorfile>
      <anchor>aa8f55d61e23edd951db0435f7c49f692</anchor>
      <arglist>(sim_solver *solver, sim_in *in, sim_out *out)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>blasfeo_dmat_args</name>
    <filename>df/dc7/structblasfeo__dmat__args.html</filename>
    <member kind="variable">
      <type>struct blasfeo_dmat *</type>
      <name>A</name>
      <anchorfile>df/dc7/structblasfeo__dmat__args.html</anchorfile>
      <anchor>a8c8aced08b9c433df69abb9061d3cd9e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ai</name>
      <anchorfile>df/dc7/structblasfeo__dmat__args.html</anchorfile>
      <anchor>a6f7ab170e2e64547b8e4b5598cd70a01</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>aj</name>
      <anchorfile>df/dc7/structblasfeo__dmat__args.html</anchorfile>
      <anchor>a9a148f7afac37d7c1912a505c2a20d0c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>blasfeo_dvec_args</name>
    <filename>d4/ddf/structblasfeo__dvec__args.html</filename>
    <member kind="variable">
      <type>struct blasfeo_dvec *</type>
      <name>x</name>
      <anchorfile>d4/ddf/structblasfeo__dvec__args.html</anchorfile>
      <anchor>a75733bc0ff85268b6a69bd87fb0076a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>xi</name>
      <anchorfile>d4/ddf/structblasfeo__dvec__args.html</anchorfile>
      <anchor>afe2c7aa282a1ba0b65fb0963cab5b862</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>colmaj_args</name>
    <filename>d5/d99/structcolmaj__args.html</filename>
    <member kind="variable">
      <type>double *</type>
      <name>A</name>
      <anchorfile>d5/d99/structcolmaj__args.html</anchorfile>
      <anchor>a2b33f8415461aefc788ce8ceb5855de1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>lda</name>
      <anchorfile>d5/d99/structcolmaj__args.html</anchorfile>
      <anchor>a8d42a4e12b27600aed0a095065e61397</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>condensing_module</name>
    <filename>d5/d7b/structcondensing__module.html</filename>
    <member kind="variable">
      <type>ocp_qp_condensing_config *</type>
      <name>config</name>
      <anchorfile>d5/d7b/structcondensing__module.html</anchorfile>
      <anchor>ac91bd1c8f7b9022a0bf7d1af134cdd6e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>dims</name>
      <anchorfile>d5/d7b/structcondensing__module.html</anchorfile>
      <anchor>a3daf105e225638bc12b50498a7eb5068</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>opts</name>
      <anchorfile>d5/d7b/structcondensing__module.html</anchorfile>
      <anchor>a7ba9bdcbbd8230bbd9dc474ab8efd16d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>mem</name>
      <anchorfile>d5/d7b/structcondensing__module.html</anchorfile>
      <anchor>a797b7f8d90404041ced6431756dbebcf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>work</name>
      <anchorfile>d5/d7b/structcondensing__module.html</anchorfile>
      <anchor>aa59a99689f6a206b0abe41142c566584</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>condensing_plan</name>
    <filename>d8/d72/structcondensing__plan.html</filename>
    <member kind="variable">
      <type>condensing_t</type>
      <name>condensing_type</name>
      <anchorfile>d8/d72/structcondensing__plan.html</anchorfile>
      <anchor>a98f3d5101ad3933b80a57107caee88e5</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dense_qp_info</name>
    <filename>d7/d84/structdense__qp__info.html</filename>
    <member kind="variable">
      <type>double</type>
      <name>solve_QP_time</name>
      <anchorfile>d7/d84/structdense__qp__info.html</anchorfile>
      <anchor>a09ad71a2f952fd3522666b3e021ad643</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>interface_time</name>
      <anchorfile>d7/d84/structdense__qp__info.html</anchorfile>
      <anchor>a8379bab37f9f162c06e36afb54cb069f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>total_time</name>
      <anchorfile>d7/d84/structdense__qp__info.html</anchorfile>
      <anchor>a92daa611308c55034919324dd5df553e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_iter</name>
      <anchorfile>d7/d84/structdense__qp__info.html</anchorfile>
      <anchor>af3d29a7ffd20409e95fcedcdc0267163</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>t_computed</name>
      <anchorfile>d7/d84/structdense__qp__info.html</anchorfile>
      <anchor>ad7e70193f4b8ffe2f5e730eaae0b1b24</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dense_qp_solver</name>
    <filename>d8/da8/structdense__qp__solver.html</filename>
    <member kind="variable">
      <type>qp_solver_config *</type>
      <name>config</name>
      <anchorfile>d8/da8/structdense__qp__solver.html</anchorfile>
      <anchor>aff7a0bdefd3f6ecb72f5f309739a05a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>dims</name>
      <anchorfile>d8/da8/structdense__qp__solver.html</anchorfile>
      <anchor>aa63b12bea81a02822c842682221ad200</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>opts</name>
      <anchorfile>d8/da8/structdense__qp__solver.html</anchorfile>
      <anchor>a2b0d3ef08a2e987c354ee64bace15baf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>mem</name>
      <anchorfile>d8/da8/structdense__qp__solver.html</anchorfile>
      <anchor>a2f4e57d0e05cba735bd5033fa4f5fc72</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>work</name>
      <anchorfile>d8/da8/structdense__qp__solver.html</anchorfile>
      <anchor>a901f9bbdc51b1d2bd4b7c5345a64d370</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>dense_qp_solver_plan</name>
    <filename>d3/df9/structdense__qp__solver__plan.html</filename>
    <member kind="variable">
      <type>dense_qp_solver_t</type>
      <name>qp_solver</name>
      <anchorfile>d3/df9/structdense__qp__solver__plan.html</anchorfile>
      <anchor>a0f2c083a2c15462eae2c55cca134f2d2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>external_function_casadi</name>
    <filename>d4/d00/structexternal__function__casadi.html</filename>
    <member kind="variable">
      <type>void(*</type>
      <name>evaluate</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>ade5f532ddc806aea69996374f1b3fcdf</anchor>
      <arglist>)(void *, ext_fun_arg_t *, void **, ext_fun_arg_t *, void **)</arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>ptr_ext_mem</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>a2ab8ec5f602fcb7f45f4993f393d7287</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>casadi_fun</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>a6797943136b545bf3abafe43521017ee</anchor>
      <arglist>)(const double **, double **, int *, double *, void *)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>casadi_work</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>a30105f1984dfbec07afe8aa44f535e7f</anchor>
      <arglist>)(int *, int *, int *, int *)</arglist>
    </member>
    <member kind="variable">
      <type>const int *(*</type>
      <name>casadi_sparsity_in</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>ac5923fd3904ba5c5d33c124a70086fb1</anchor>
      <arglist>)(int)</arglist>
    </member>
    <member kind="variable">
      <type>const int *(*</type>
      <name>casadi_sparsity_out</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>a87282c6fc8adb049b28e9fe004efcefa</anchor>
      <arglist>)(int)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>casadi_n_in</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>a9170b97bd8689a7d3ab98cb4bf35af8d</anchor>
      <arglist>)()</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>casadi_n_out</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>a2afec380039c854870fb8425bfc2b3bc</anchor>
      <arglist>)()</arglist>
    </member>
    <member kind="variable">
      <type>double **</type>
      <name>args</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>a92a6553705e2d7881892f79cd0210a21</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double **</type>
      <name>res</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>a06b6460bd9ecfadba1bfe7261dfd253c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>w</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>ae9bf29d6890c9b131d0b0ec96f410942</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>iw</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>a0bc3c11c37acfad25c782c04b0acd6c7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>args_size</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>a96bf63d74aa0ea1973230a42072293b5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>res_size</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>ac756283aca0aa1ef5fcc811e09545b6a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>args_num</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>a21107be798588d4ad93664085ae09192</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>args_size_tot</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>a6accc3cc231122f8d4c8bbba5a143690</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>res_num</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>ad8aae9f8ee09a07a3daf457dfcf59b5d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>res_size_tot</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>afb8dd9ce7285cc28b98a5dd02aaed187</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>in_num</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>a7985d88140e62b0af3025385b5e0cf96</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>out_num</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>ad7fb0cb5e9ef1777c92223f9865179d2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>iw_size</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>a81430ef4c8904b6024ba09083c8249fb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>w_size</name>
      <anchorfile>d4/d00/structexternal__function__casadi.html</anchorfile>
      <anchor>ada522f8ec9f506a8e492792eb682caf9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>external_function_generic</name>
    <filename>d0/da6/structexternal__function__generic.html</filename>
    <member kind="variable">
      <type>void(*</type>
      <name>evaluate</name>
      <anchorfile>d0/da6/structexternal__function__generic.html</anchorfile>
      <anchor>aa9ba52cbe303decbb17c1ecbab1124f6</anchor>
      <arglist>)(void *, ext_fun_arg_t *, void **, ext_fun_arg_t *, void **)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>external_function_param_casadi</name>
    <filename>d9/d4d/structexternal__function__param__casadi.html</filename>
    <member kind="variable">
      <type>void(*</type>
      <name>evaluate</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>a5d1b54592f635b4eb2e171d105372d48</anchor>
      <arglist>)(void *, ext_fun_arg_t *, void **, ext_fun_arg_t *, void **)</arglist>
    </member>
    <member kind="variable">
      <type>void(*</type>
      <name>set_param</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>a56db5b2804ebb289cd3e8f892c042b21</anchor>
      <arglist>)(void *, double *)</arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>ptr_ext_mem</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>a0ac139eb850ed7bf65f412cfaa28e3ce</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>casadi_fun</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>a3a1e26bb5fc120d9c0d12e2875dd3fa5</anchor>
      <arglist>)(const double **, double **, int *, double *, void *)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>casadi_work</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>afdbde95e03a312407d0734ec04d1e3a5</anchor>
      <arglist>)(int *, int *, int *, int *)</arglist>
    </member>
    <member kind="variable">
      <type>const int *(*</type>
      <name>casadi_sparsity_in</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>a899d308c3219708afd6615d1017e2d0c</anchor>
      <arglist>)(int)</arglist>
    </member>
    <member kind="variable">
      <type>const int *(*</type>
      <name>casadi_sparsity_out</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>aed65e7a569584df40abb5c4cb4211467</anchor>
      <arglist>)(int)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>casadi_n_in</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>a356dc4c83f945d7e68e4a9120e635d80</anchor>
      <arglist>)()</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>casadi_n_out</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>aadb1f42b0b9cec3724e7839fcd9d892d</anchor>
      <arglist>)()</arglist>
    </member>
    <member kind="variable">
      <type>double **</type>
      <name>args</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>abc7e577e6d376f38ceaf4c728e10ef5d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double **</type>
      <name>res</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>a554d3260c5d60f79a7a99aad9d1768e1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>w</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>a8c0710d01954c9f4c000136d2132514b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>p</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>a972e5dd6b476287a02f3dbdec65453c8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>iw</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>a2f2c8f9ea7f876284d1dac627a6824f6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>args_size</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>abaf580e900737a017932f24776060e6a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>res_size</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>a17223ce92404ab32e3f67391af718ba4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>args_num</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>a2005173bd53b4a4f6a9b6f6db0eda208</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>args_size_tot</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>a21b5304313d2a6ee45f8bae20f67f2dd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>res_num</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>aa1ada92ab6827d214e1714535602dff8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>res_size_tot</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>ab8afc8074dcd882f57dfeb5565f20f67</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>in_num</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>ad6fb3a4989e74c3349efd8846d5efade</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>out_num</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>a91d1804fc202dffa4073846c4a4d7fca</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>iw_size</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>ad42a2e1a5d2d6745658e0195b6740b3d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>w_size</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>ad8e3c6ff86968c8e4366964558cb1b25</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>np</name>
      <anchorfile>d9/d4d/structexternal__function__param__casadi.html</anchorfile>
      <anchor>adf59d31d9d251011803f71c25352e085</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>ocp_nlp_plan</name>
    <filename>db/d0e/structocp__nlp__plan.html</filename>
    <member kind="variable">
      <type>ocp_qp_solver_plan</type>
      <name>ocp_qp_solver_plan</name>
      <anchorfile>db/d0e/structocp__nlp__plan.html</anchorfile>
      <anchor>a503dffd0eecf0ea812c2e6bb1fe9617b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>sim_solver_plan *</type>
      <name>sim_solver_plan</name>
      <anchorfile>db/d0e/structocp__nlp__plan.html</anchorfile>
      <anchor>a09d2886157409a50f2211a569047c231</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>ocp_nlp_solver_t</type>
      <name>nlp_solver</name>
      <anchorfile>db/d0e/structocp__nlp__plan.html</anchorfile>
      <anchor>a0e775617625534cda4ca933ae673e4cc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>ocp_nlp_reg_t</type>
      <name>regularization</name>
      <anchorfile>db/d0e/structocp__nlp__plan.html</anchorfile>
      <anchor>aaa527d9902ede15e1c48230ed3003653</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>ocp_nlp_cost_t *</type>
      <name>nlp_cost</name>
      <anchorfile>db/d0e/structocp__nlp__plan.html</anchorfile>
      <anchor>a2e2811f79052dab8b5265538eb8bf711</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>ocp_nlp_dynamics_t *</type>
      <name>nlp_dynamics</name>
      <anchorfile>db/d0e/structocp__nlp__plan.html</anchorfile>
      <anchor>a8d1dffdd081b06252cbb0be8c349f5b0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>ocp_nlp_constraints_t *</type>
      <name>nlp_constraints</name>
      <anchorfile>db/d0e/structocp__nlp__plan.html</anchorfile>
      <anchor>a081148f501f71f0aa9823d3ea87ba73f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>N</name>
      <anchorfile>db/d0e/structocp__nlp__plan.html</anchorfile>
      <anchor>a724e8ea8260a803bb4767c3d00a4d7b1</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>ocp_nlp_solver</name>
    <filename>d2/de0/structocp__nlp__solver.html</filename>
    <member kind="variable">
      <type>ocp_nlp_config *</type>
      <name>config</name>
      <anchorfile>d2/de0/structocp__nlp__solver.html</anchorfile>
      <anchor>a9c482395f5105894dc788f8da5313e90</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>dims</name>
      <anchorfile>d2/de0/structocp__nlp__solver.html</anchorfile>
      <anchor>aec54fc5ba2d9ed0aacb0c74bea2b3479</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>opts</name>
      <anchorfile>d2/de0/structocp__nlp__solver.html</anchorfile>
      <anchor>ae4e0fc0e95b16ab0157d51516677e964</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>mem</name>
      <anchorfile>d2/de0/structocp__nlp__solver.html</anchorfile>
      <anchor>aca36b3846ee37a5d0180358febb62ad0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>work</name>
      <anchorfile>d2/de0/structocp__nlp__solver.html</anchorfile>
      <anchor>a3cfb1a97a99e29531915a474079cba78</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>ocp_qp_condensing_config</name>
    <filename>da/d5e/structocp__qp__condensing__config.html</filename>
    <member kind="variable">
      <type>int(*</type>
      <name>condensing</name>
      <anchorfile>da/d5e/structocp__qp__condensing__config.html</anchorfile>
      <anchor>af0100f5710f15be3d0d359749c2aae8b</anchor>
      <arglist>)(void *qp_in, void *qp_out, void *opts, void *mem, void *work)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>expansion</name>
      <anchorfile>da/d5e/structocp__qp__condensing__config.html</anchorfile>
      <anchor>af43bd7eddc43bc006af61fae49b84fd6</anchor>
      <arglist>)(void *qp_in, void *qp_out, void *opts, void *mem, void *work)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>opts_calculate_size</name>
      <anchorfile>da/d5e/structocp__qp__condensing__config.html</anchorfile>
      <anchor>a4e35de85149b1b76aa162d36b5060145</anchor>
      <arglist>)(ocp_qp_dims *dims)</arglist>
    </member>
    <member kind="variable">
      <type>void *(*</type>
      <name>opts_assign</name>
      <anchorfile>da/d5e/structocp__qp__condensing__config.html</anchorfile>
      <anchor>aef92853425c2dc55b3345f556f9c3f18</anchor>
      <arglist>)(ocp_qp_dims *dims, void *raw_memory)</arglist>
    </member>
    <member kind="variable">
      <type>void(*</type>
      <name>opts_initialize_default</name>
      <anchorfile>da/d5e/structocp__qp__condensing__config.html</anchorfile>
      <anchor>ad3bb8514481854ad9075cfae3c47512c</anchor>
      <arglist>)(ocp_qp_dims *dims, void *opts)</arglist>
    </member>
    <member kind="variable">
      <type>void(*</type>
      <name>opts_update</name>
      <anchorfile>da/d5e/structocp__qp__condensing__config.html</anchorfile>
      <anchor>a11057f8650e3cd4d6ea5e083d11f4047</anchor>
      <arglist>)(ocp_qp_dims *dims, void *opts)</arglist>
    </member>
    <member kind="variable">
      <type>void(*</type>
      <name>opts_set</name>
      <anchorfile>da/d5e/structocp__qp__condensing__config.html</anchorfile>
      <anchor>a925f8e35b32db5cf08df0bb3f11d26a7</anchor>
      <arglist>)(void *config_, void *opts_, const char *field, const void *value)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>memory_calculate_size</name>
      <anchorfile>da/d5e/structocp__qp__condensing__config.html</anchorfile>
      <anchor>a46cd19de78d7b7d2c0b8dca248e0340a</anchor>
      <arglist>)(ocp_qp_dims *dims, void *opts)</arglist>
    </member>
    <member kind="variable">
      <type>void *(*</type>
      <name>memory_assign</name>
      <anchorfile>da/d5e/structocp__qp__condensing__config.html</anchorfile>
      <anchor>a462844cf24610295468cdf226e3fa176</anchor>
      <arglist>)(ocp_qp_dims *dims, void *opts, void *raw_memory)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>workspace_calculate_size</name>
      <anchorfile>da/d5e/structocp__qp__condensing__config.html</anchorfile>
      <anchor>ac343ca679d074b4bfd5cd2c87edd09b0</anchor>
      <arglist>)(ocp_qp_dims *dims, void *opts)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>ocp_qp_info</name>
    <filename>d1/d83/structocp__qp__info.html</filename>
    <member kind="variable">
      <type>double</type>
      <name>solve_QP_time</name>
      <anchorfile>d1/d83/structocp__qp__info.html</anchorfile>
      <anchor>a2507c28d8e718153c4292ada11dcd0b6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>condensing_time</name>
      <anchorfile>d1/d83/structocp__qp__info.html</anchorfile>
      <anchor>a05ee07c650a53d97c7a594fa4cdfabc0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>interface_time</name>
      <anchorfile>d1/d83/structocp__qp__info.html</anchorfile>
      <anchor>afe86668904c93655974d4cb6c1fa8478</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>total_time</name>
      <anchorfile>d1/d83/structocp__qp__info.html</anchorfile>
      <anchor>a93a86098f8161b085fdda7d4313dda8c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_iter</name>
      <anchorfile>d1/d83/structocp__qp__info.html</anchorfile>
      <anchor>a35f1847736d313a27381fe76f760fdbe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>t_computed</name>
      <anchorfile>d1/d83/structocp__qp__info.html</anchorfile>
      <anchor>afcdcab629872e40910bc35071dd4d1b7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>ocp_qp_solver</name>
    <filename>d2/d16/structocp__qp__solver.html</filename>
    <member kind="variable">
      <type>ocp_qp_xcond_solver_config *</type>
      <name>config</name>
      <anchorfile>d2/d16/structocp__qp__solver.html</anchorfile>
      <anchor>a21e7474c8d5334d45c471c59ea43d1e8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>dims</name>
      <anchorfile>d2/d16/structocp__qp__solver.html</anchorfile>
      <anchor>a7afdfa5d7f15817202e1a284df85f59d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>opts</name>
      <anchorfile>d2/d16/structocp__qp__solver.html</anchorfile>
      <anchor>a83c4efb7a0a27e1ba2565df6bbae6d13</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>mem</name>
      <anchorfile>d2/d16/structocp__qp__solver.html</anchorfile>
      <anchor>ac97674ad586fff6794030c9901bd2481</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>work</name>
      <anchorfile>d2/d16/structocp__qp__solver.html</anchorfile>
      <anchor>aaa461d90307cf121a632c45c15a08faa</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>ocp_qp_solver_plan</name>
    <filename>d1/d17/structocp__qp__solver__plan.html</filename>
    <member kind="variable">
      <type>ocp_qp_solver_t</type>
      <name>qp_solver</name>
      <anchorfile>d1/d17/structocp__qp__solver__plan.html</anchorfile>
      <anchor>affb22bb0053b8a7fcbfb7a458522ed3d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>ocp_qp_xcond_solver_config</name>
    <filename>d7/d86/structocp__qp__xcond__solver__config.html</filename>
    <member kind="variable">
      <type>void(*</type>
      <name>dims_set</name>
      <anchorfile>d7/d86/structocp__qp__xcond__solver__config.html</anchorfile>
      <anchor>a77d424c6b9f049a97cd4c45bf3ed9850</anchor>
      <arglist>)(void *config_, void *dims_, int stage, const char *field, const int *value)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>evaluate</name>
      <anchorfile>d7/d86/structocp__qp__xcond__solver__config.html</anchorfile>
      <anchor>aafaeb10b583afd85848f58da86641550</anchor>
      <arglist>)(void *config, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts, void *mem, void *work)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>opts_calculate_size</name>
      <anchorfile>d7/d86/structocp__qp__xcond__solver__config.html</anchorfile>
      <anchor>ac2f8c70d6cd40a16bd95f57b11fcd0c1</anchor>
      <arglist>)(void *config, ocp_qp_dims *dims)</arglist>
    </member>
    <member kind="variable">
      <type>void *(*</type>
      <name>opts_assign</name>
      <anchorfile>d7/d86/structocp__qp__xcond__solver__config.html</anchorfile>
      <anchor>a422fbb1f04e643d1cbd23f3d7e620b7b</anchor>
      <arglist>)(void *config, ocp_qp_dims *dims, void *raw_memory)</arglist>
    </member>
    <member kind="variable">
      <type>void(*</type>
      <name>opts_initialize_default</name>
      <anchorfile>d7/d86/structocp__qp__xcond__solver__config.html</anchorfile>
      <anchor>ae74fde235301417ec76cd79bbd0a9991</anchor>
      <arglist>)(void *config, ocp_qp_dims *dims, void *opts)</arglist>
    </member>
    <member kind="variable">
      <type>void(*</type>
      <name>opts_update</name>
      <anchorfile>d7/d86/structocp__qp__xcond__solver__config.html</anchorfile>
      <anchor>a607a93747ae381485cefe1657bdeb6d7</anchor>
      <arglist>)(void *config, ocp_qp_dims *dims, void *opts)</arglist>
    </member>
    <member kind="variable">
      <type>void(*</type>
      <name>opts_set</name>
      <anchorfile>d7/d86/structocp__qp__xcond__solver__config.html</anchorfile>
      <anchor>ad6782a69356e7fb583f0c1df0b2f2e17</anchor>
      <arglist>)(void *config_, void *opts_, const char *field, const void *value)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>memory_calculate_size</name>
      <anchorfile>d7/d86/structocp__qp__xcond__solver__config.html</anchorfile>
      <anchor>aae4df50d30279ecdbe6a4e98aa934206</anchor>
      <arglist>)(void *config, ocp_qp_dims *dims, void *opts)</arglist>
    </member>
    <member kind="variable">
      <type>void *(*</type>
      <name>memory_assign</name>
      <anchorfile>d7/d86/structocp__qp__xcond__solver__config.html</anchorfile>
      <anchor>a1f3f4d781e480ad24ed9372b9a1e688c</anchor>
      <arglist>)(void *config, ocp_qp_dims *dims, void *opts, void *raw_memory)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>workspace_calculate_size</name>
      <anchorfile>d7/d86/structocp__qp__xcond__solver__config.html</anchorfile>
      <anchor>a7687042a5a720a97a99836af971497c2</anchor>
      <arglist>)(void *config, ocp_qp_dims *dims, void *opts)</arglist>
    </member>
    <member kind="variable">
      <type>qp_solver_config *</type>
      <name>qp_solver</name>
      <anchorfile>d7/d86/structocp__qp__xcond__solver__config.html</anchorfile>
      <anchor>a5c6090a7676a76f925745281a0e99305</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>qp_solver_config</name>
    <filename>d3/db2/structqp__solver__config.html</filename>
    <member kind="variable">
      <type>void(*</type>
      <name>dims_set</name>
      <anchorfile>d3/db2/structqp__solver__config.html</anchorfile>
      <anchor>aa6067a91465ce85442cc2f2e62e5645a</anchor>
      <arglist>)(void *config_, void *dims_, int stage, const char *field, const int *value)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>evaluate</name>
      <anchorfile>d3/db2/structqp__solver__config.html</anchorfile>
      <anchor>add435247e88058363505a693bec0e991</anchor>
      <arglist>)(void *config, void *qp_in, void *qp_out, void *opts, void *mem, void *work)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>opts_calculate_size</name>
      <anchorfile>d3/db2/structqp__solver__config.html</anchorfile>
      <anchor>a0147c5417a9f3fad015f7266f39a2b6a</anchor>
      <arglist>)(void *config, void *dims)</arglist>
    </member>
    <member kind="variable">
      <type>void *(*</type>
      <name>opts_assign</name>
      <anchorfile>d3/db2/structqp__solver__config.html</anchorfile>
      <anchor>a8f590dc1a06d25c900ba1d623f33c1b6</anchor>
      <arglist>)(void *config, void *dims, void *raw_memory)</arglist>
    </member>
    <member kind="variable">
      <type>void(*</type>
      <name>opts_initialize_default</name>
      <anchorfile>d3/db2/structqp__solver__config.html</anchorfile>
      <anchor>a89c4e1d2c86ead2aacddf68b2e4ba43f</anchor>
      <arglist>)(void *config, void *dims, void *opts)</arglist>
    </member>
    <member kind="variable">
      <type>void(*</type>
      <name>opts_update</name>
      <anchorfile>d3/db2/structqp__solver__config.html</anchorfile>
      <anchor>aa3a666c91c91588e2001d42010f68fe8</anchor>
      <arglist>)(void *config, void *dims, void *opts)</arglist>
    </member>
    <member kind="variable">
      <type>void(*</type>
      <name>opts_set</name>
      <anchorfile>d3/db2/structqp__solver__config.html</anchorfile>
      <anchor>a2eb192932af0a52dfee9e56ec3887fc9</anchor>
      <arglist>)(void *config_, void *opts_, const char *field, const void *value)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>memory_calculate_size</name>
      <anchorfile>d3/db2/structqp__solver__config.html</anchorfile>
      <anchor>ad8ed4ae79f6e5b258520d16234a52561</anchor>
      <arglist>)(void *config, void *dims, void *opts)</arglist>
    </member>
    <member kind="variable">
      <type>void *(*</type>
      <name>memory_assign</name>
      <anchorfile>d3/db2/structqp__solver__config.html</anchorfile>
      <anchor>a8d19640e720a73ae61ff6331a561cc69</anchor>
      <arglist>)(void *config, void *dims, void *opts, void *raw_memory)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>workspace_calculate_size</name>
      <anchorfile>d3/db2/structqp__solver__config.html</anchorfile>
      <anchor>a7de2a9e155ac1942c228873b367aca5f</anchor>
      <arglist>)(void *config, void *dims, void *opts)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>sim_config</name>
    <filename>d0/d7e/structsim__config.html</filename>
    <member kind="variable">
      <type>int(*</type>
      <name>evaluate</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>af814d56b1ea67adbe7d99cdc327c2e10</anchor>
      <arglist>)(void *config_, sim_in *in, sim_out *out, void *opts, void *mem, void *work)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>precompute</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>afbbb6994ee76ea6288ab062ae578ffd0</anchor>
      <arglist>)(void *config_, sim_in *in, sim_out *out, void *opts, void *mem, void *work)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>opts_calculate_size</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>afbe7b347c6381f7bd1d9daf5d080e5f8</anchor>
      <arglist>)(void *config_, void *dims)</arglist>
    </member>
    <member kind="variable">
      <type>void *(*</type>
      <name>opts_assign</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>a37fabed0857b70473eb5746934ab5f62</anchor>
      <arglist>)(void *config_, void *dims, void *raw_memory)</arglist>
    </member>
    <member kind="variable">
      <type>void(*</type>
      <name>opts_initialize_default</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>afdc588ed7bfa4ad08dfbc9aa986c22bb</anchor>
      <arglist>)(void *config_, void *dims, void *opts)</arglist>
    </member>
    <member kind="variable">
      <type>void(*</type>
      <name>opts_update</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>a217c1217ba71c79d8b6698d12f0f05c2</anchor>
      <arglist>)(void *config_, void *dims, void *opts)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>opts_set</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>ad7f296cf6ec0e5b5fed431a03760ad44</anchor>
      <arglist>)(void *config_, void *opts_, const char *field, void *value)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>memory_calculate_size</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>a6a8ecf2263944eae6fe323ab33d7f72a</anchor>
      <arglist>)(void *config, void *dims, void *opts)</arglist>
    </member>
    <member kind="variable">
      <type>void *(*</type>
      <name>memory_assign</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>af4242ffc3fef02dff310d82ad5c15029</anchor>
      <arglist>)(void *config, void *dims, void *opts, void *raw_memory)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>workspace_calculate_size</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>a8011aaa8e1e97edfaa66ce08462437d5</anchor>
      <arglist>)(void *config, void *dims, void *opts)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>model_calculate_size</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>ae443210a9a0eb03c9f6ba83222acbfec</anchor>
      <arglist>)(void *config, void *dims)</arglist>
    </member>
    <member kind="variable">
      <type>void *(*</type>
      <name>model_assign</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>a6e5e0f4a830ef1a842af6c81aa55af57</anchor>
      <arglist>)(void *config, void *dims, void *raw_memory)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>model_set</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>a462c55f3f57ad8ccc9b2170e93ce06fb</anchor>
      <arglist>)(void *model, const char *field, void *value)</arglist>
    </member>
    <member kind="variable">
      <type>void(*</type>
      <name>config_initialize_default</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>ae522c287bf2f049513270e74c2e07423</anchor>
      <arglist>)(void *config)</arglist>
    </member>
    <member kind="variable">
      <type>int(*</type>
      <name>dims_calculate_size</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>a0b72c70ca223d62ff6b19fe4a718c8b4</anchor>
      <arglist>)()</arglist>
    </member>
    <member kind="variable">
      <type>void *(*</type>
      <name>dims_assign</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>af3f06fb1150f6ace8b2b7f8271da2de2</anchor>
      <arglist>)(void *config, void *raw_memory)</arglist>
    </member>
    <member kind="variable">
      <type>void(*</type>
      <name>dims_set</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>a86bb381023726322f0353757adbcbb97</anchor>
      <arglist>)(void *config, void *dims, const char *field, const int *value)</arglist>
    </member>
    <member kind="variable">
      <type>void(*</type>
      <name>dims_get</name>
      <anchorfile>d0/d7e/structsim__config.html</anchorfile>
      <anchor>a89eb9499fbf2cd5bb45f1fed16ed2a3b</anchor>
      <arglist>)(void *config, void *dims, const char *field, int *value)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>sim_in</name>
    <filename>d0/ddd/structsim__in.html</filename>
    <member kind="variable">
      <type>void *</type>
      <name>dims</name>
      <anchorfile>d0/ddd/structsim__in.html</anchorfile>
      <anchor>a085794c1443de1af19f007559a3e87a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>x</name>
      <anchorfile>d0/ddd/structsim__in.html</anchorfile>
      <anchor>ac23f5c3f1115c1173029b1d166b681d2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>u</name>
      <anchorfile>d0/ddd/structsim__in.html</anchorfile>
      <anchor>a0025c3cba13b54768dad01ae3677c01e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>xdot</name>
      <anchorfile>d0/ddd/structsim__in.html</anchorfile>
      <anchor>ab7a0d01aeb5251c01270338880950ea1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>z</name>
      <anchorfile>d0/ddd/structsim__in.html</anchorfile>
      <anchor>a4c6726311bbb0152d0c3fd53a0dd9475</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>S_forw</name>
      <anchorfile>d0/ddd/structsim__in.html</anchorfile>
      <anchor>ac5bb55e037a908797dabd07ae17d7288</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>S_adj</name>
      <anchorfile>d0/ddd/structsim__in.html</anchorfile>
      <anchor>ac87faee3738dd7bd49e39094a8652565</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>model</name>
      <anchorfile>d0/ddd/structsim__in.html</anchorfile>
      <anchor>a424c9c66ed9b5b2e21a5c7aac20b13db</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>T</name>
      <anchorfile>d0/ddd/structsim__in.html</anchorfile>
      <anchor>a54f85591d8934415d5f5d66a697d5a29</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>sim_info</name>
    <filename>db/d02/structsim__info.html</filename>
    <member kind="variable">
      <type>double</type>
      <name>CPUtime</name>
      <anchorfile>db/d02/structsim__info.html</anchorfile>
      <anchor>a2967d18408908859de93022c447d7d7e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>LAtime</name>
      <anchorfile>db/d02/structsim__info.html</anchorfile>
      <anchor>a60d536daa57dab1a77adcfd289a951c7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ADtime</name>
      <anchorfile>db/d02/structsim__info.html</anchorfile>
      <anchor>ab2def1f477def553a0ecee2e20fd6de3</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>sim_opts</name>
    <filename>dc/dd2/structsim__opts.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>ns</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>a5ac298a18e51680d7218df919b919a5e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_steps</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>a3bf07bb203769c48af0e15804fa2ffff</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_forw_sens</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>a9aadd5674e2a93a0f54898673a8ca29c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>tableau_size</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>ae6bdc8af312aa40190082dbc2e32c03c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>A_mat</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>a1c7c96075c738ff2f211230727ebb082</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>c_vec</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>a4944cbfc40941f39fa990af112bcd346</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>b_vec</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>a621b740e9db21184f9f365ba11a759d3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>sens_forw</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>aa137aaf37374e4d59a3b294d09a0cabd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>sens_adj</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>a10f3023d1bbeec53f9aeebc15707e574</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>sens_hess</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>afb8d6e2179d8e3ed2192bd6289066bc9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>output_z</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>aec2e51f80897da2f56b0218c74db5852</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>sens_algebraic</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>a695d585a493128157197a5fe043fb622</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>newton_iter</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>a378d17cc8d33f36cd6869c450e6d9797</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>jac_reuse</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>a4d48c780bfed35c58d1d79e0c2dc022a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Newton_scheme *</type>
      <name>scheme</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>a63627e83d2b186b664accb21d8e92474</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>work</name>
      <anchorfile>dc/dd2/structsim__opts.html</anchorfile>
      <anchor>a3d297030ff7528cafbbb46e91d012235</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>sim_out</name>
    <filename>db/d9b/structsim__out.html</filename>
    <member kind="variable">
      <type>double *</type>
      <name>xn</name>
      <anchorfile>db/d9b/structsim__out.html</anchorfile>
      <anchor>a3c24d3363a34bf9be4dc8fc8b7531f22</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>S_forw</name>
      <anchorfile>db/d9b/structsim__out.html</anchorfile>
      <anchor>ab4953f90222c679ce79f4f0ba9599294</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>S_adj</name>
      <anchorfile>db/d9b/structsim__out.html</anchorfile>
      <anchor>aa9eac3f3cced012c9533bc34a5df2edc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>S_hess</name>
      <anchorfile>db/d9b/structsim__out.html</anchorfile>
      <anchor>a42cbaff33ce2ab39afcf804f1440477a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>zn</name>
      <anchorfile>db/d9b/structsim__out.html</anchorfile>
      <anchor>a162d8eec101d71f38cec0068a2b6e89c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>S_algebraic</name>
      <anchorfile>db/d9b/structsim__out.html</anchorfile>
      <anchor>ab1a9fcae6fb5090386e9520080cb3093</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>grad</name>
      <anchorfile>db/d9b/structsim__out.html</anchorfile>
      <anchor>ac71c395630303f1b9747ebf5517ad004</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>sim_info *</type>
      <name>info</name>
      <anchorfile>db/d9b/structsim__out.html</anchorfile>
      <anchor>ab6b88c52b7487ea35b2c27f64b1e670a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>sim_solver</name>
    <filename>df/d96/structsim__solver.html</filename>
    <member kind="variable">
      <type>sim_config *</type>
      <name>config</name>
      <anchorfile>df/d96/structsim__solver.html</anchorfile>
      <anchor>a0051c658bc3b2685507a7b49aad063b1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>dims</name>
      <anchorfile>df/d96/structsim__solver.html</anchorfile>
      <anchor>a5a2641cc852eef28a6497397343ad501</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>opts</name>
      <anchorfile>df/d96/structsim__solver.html</anchorfile>
      <anchor>a7e93c7e0da2dd312859fb117f708cb5c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>mem</name>
      <anchorfile>df/d96/structsim__solver.html</anchorfile>
      <anchor>a2cff7053e8b5a9f4717326ddbab63b2b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>work</name>
      <anchorfile>df/d96/structsim__solver.html</anchorfile>
      <anchor>a521c0d25883b2136bf608a9f65c6c809</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>sim_solver_plan</name>
    <filename>d7/dd4/structsim__solver__plan.html</filename>
    <member kind="variable">
      <type>sim_solver_t</type>
      <name>sim_solver</name>
      <anchorfile>d7/dd4/structsim__solver__plan.html</anchorfile>
      <anchor>a37fd040f4cab2ee9ab8b75cf5dc99104</anchor>
      <arglist></arglist>
    </member>
  </compound>
</tagfile>
