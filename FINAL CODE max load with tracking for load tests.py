# Full simulation set up

# Import standard libraries

import numpy as np
from scipy.stats import t
from scipy.stats import norm # add brittle check for redistribution and track ductility of each potential hinge location, check how neutral axis is being determined at different stages
from scipy.integrate import simpson
import matplotlib.pyplot as plt
import scipy.stats as stats

# Installed FEA libraries

from Pynite import FEModel3D 
from Pynite.Visualization import Renderer # only required for visual verification

# FEA function

def FEA(E_c_dict, I_dict, k_dict, scaled_node_load_dict):
    
    # Create FE model

    model_1 = FEModel3D()

    # Generate model nodes

    for i in range (num_grid):
        
        if i % 2 != 0  and num_points <= i < num_points * 2: 
        
           x_coord = (i - num_points) * rib_spacing / 2
           y_coord = 0
           z_coord = 0
        
           node_name = f'node_{i}'
        
           model_1.add_node(node_name, x_coord, y_coord, z_coord)
   
    # Create beam sections

    b_model = 100 # arbitrary values
    h_model = 100 

    A_model = b_model * h_model 
    Iy_model = (b_model * h_model ** 3) / 12
    Iz_model = (h_model * b_model ** 3) / 12
    J_model = ((b_model * h_model) / 12 * (b_model**2 + h_model**2))

    section_name = 'section'

    model_1.add_section(section_name, A_model, Iy_model, Iz_model, J_model)
        
    # Create beam materials

    for i in range(num_grid):
        
        if i % 2 == 0  and num_points <= i < num_points * 2:
        
           E_i = reduction_factor * E_c_dict[i] * I_dict[i] / Iy_model
           nu_i = poissons_ratio
           G_i = E_i / (2 * (1 + nu_i))
           rho_i = 1 # arbitrary as self weight is calculated seperately
        
           material_name = f'material_{i}'
       
           model_1.add_material(material_name, E_i, G_i, nu_i, rho_i)
        
    for i in range(num_grid):
       
        if i % 2 == 0  and num_points <= i < num_points * 2 - 1: 
           
           start_node = f'node_{i-1}'
           end_node = f'node_{i+1}'
        
           member_name = f'member_{i}'
           material_name = f'material_{i}'
           section_name = 'section'
        
           model_1.add_member(member_name, start_node, end_node, material_name, section_name)
  
    # Apply boundary conditions using spring supports

    for i in range(num_grid):
        
        if i % 2 != 0  and num_points <= i < num_points * 2: 
        
           node_name = f'node_{i}'
        
           k_support = k_dict[i]
           
          
      
           ##model_1.def_support(node_name, True, True, True, False, False, False)
           
           model_1.def_support_spring(node_name, 'DY', k_support)
      
    # Restrict displacements and rotations for edge nodes
        
    pin_node_1 = f'node_{num_points}'
    pin_node_2 = f'node_{num_grid - num_points - 1}'

    model_1.def_support(pin_node_1, True, False, True, True, True, False)
    model_1.def_support(pin_node_2, True, False, True, True, True, False)

    # Apply node load to each support

    for i in range(num_grid):
         
         if i % 2 != 0  and num_points <= i < num_points * 2:       
          
           node_name = f'node_{i}'
        
           node_load = -1 * scaled_node_load_dict[i]
           
          
        
           model_1.add_node_load(node_name, 'FY', node_load, case = 'F')
        
    # Create load combination

    model_1.add_load_combo('combo', {'F':1.0})
    
    return model_1

def rib_output(k_dict, rib_ductile_fail_dict, rib_brittle_fail_dict, rib_status_dict, model_1, strain_c_ult_dict, System_ULS_fail_counter, x_dict, M_capacity_dict, EI_eq_dict, d_dict, strain_s_ult_dict, plastic_hinge_end_1_dict, plastic_hinge_end_2_dict, plastic_hinge_mid_dict, v_capacity_dict, b_dict, z_dict, h_f_dict, tau_capacity_dict, f_c_dict, M_capacity_ult_dict, deflection_max_dict, sec_brittle_dict, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter):
    
    # Extract nodal forces
    
    for i in range(num_grid):
         
        if i % 2 != 0  and num_points <= i < num_points * 2:
    
            node_name = f'node_{i}'
         
            F_float = model_1.nodes[node_name].RxnFY
       
            F = abs(F_float['combo']) / span # force on each rib as a udl
            
            M_demand_end = ((F * (span ** 2)) / 12) # calculate bending moment
            
            M_demand_end_dict[i] = M_demand_end
            
            M_demand_mid = ((F * (span ** 2)) / 24) # calculate bending moment

            
            M_demand_mid_dict[i] = M_demand_mid
         
            V_demand = (F * span) / 2 # calculate shear force
            
            V_demand_dict[i] = V_demand
            
            deflection_float = model_1.nodes[node_name].DY
            
            deflection = abs(deflection_float['combo'])
            
            deflection_dict[i] = deflection
         
    for i in range(num_grid):
        
        if i % 2 != 0  and num_points <= i < num_points * 2:
            
            M_demand_end = M_demand_end_dict[i]
         
            M_demand_mid = M_demand_mid_dict[i]
         
            M_capacity_end_1 = M_capacity_dict[i - num_points]
    
            M_capacity_end_2 = M_capacity_dict[i + num_points]
            
            M_capacity_mid = M_capacity_dict[i]
          
            M_capacity_ult_mid = M_capacity_ult_dict[i]
            
            M_capacity_ult_end_1 = M_capacity_ult_dict[i - num_points]
            
            M_capacity_ult_end_2 = M_capacity_ult_dict[i + num_points]
          
            EI = EI_eq_dict[i]
          
            M_redis_1 = 0
            M_redis_2 = 0
            M_redis_mid = 0
            
            # Check for vertical shear
        
            z = z_dict[i]
           
            v_capacity_end_1 = v_capacity_dict[i - num_points]
            
            v_capacity_end_2 = v_capacity_dict[i + num_points]
        
            v_demand = V_demand / (0.9 * rib_spacing * z) # approximation
        
            # Check end 1
            
            if v_demand > v_capacity_end_1:
                
                #print(f"Rib_{((i-num_points)/2)+1}: Vertical shear failure at end 1")
                
                rib_brittle_fail_dict[i] = True
                
                failed_rib_info.append((i, original_EI[i]))
                
                Vshear_counter +=1
                
                return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
                
            
                
            # Check end 2
            
            if v_demand > v_capacity_end_2:
                
               # print(f"Rib_{((i-num_points)/2)+1}: Vertical shear failure at end 2")
                
                Vshear_counter +=1
                
                rib_brittle_fail_dict[i] = True
                
                failed_rib_info.append((i, original_EI[i]))
                
                return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
                
             
                
            # Check longitudinal shear at mid span Eurocode 2
            
            tau_capacity = tau_capacity_dict[i]
            
            h_f = h_f_dict[i]
            
            f_c = f_c_dict[i] * conf_factor
            
            f_y = S_nom # simplicity for now
            
            delta_x = l_0_factor * 0.5 * span / 2
            
            delta_F = M_demand_mid / z
            
            tau_ed = delta_F / (h_f * delta_x)
            
            tau_rd = 0.160 * f_c * (1 - f_c / 250)
          
            if tau_rd > tau_ed:
                
                cot_theta_f = 2.0
                
            if tau_rd <= tau_ed:
                
                arg = np.clip(tau_ed / (0.2 * f_c * (1 - f_c / 250)), -1, 1)
                
                theta_f = 0.5 * np.arcsin(arg)
              
                cot_theta_f = 1 / np.tan(theta_f)
          
            tau_demand = (tau_ed * h_f) / (f_y * cot_theta_f)
            
            if tau_demand > tau_capacity:
                
               # print(f"Rib_{((i-num_points)/2)+1}: Longitudinal shear failure")
                
                rib_brittle_fail_dict[i] = True
                
                failed_rib_info.append((i, original_EI[i]))
                
                Lshear_counter += 1
                
                return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
                
          
            
            # Check brittle sections for ultimate failure (they can have moment distributed to them but cannot redistribute moment themselves)
            
            # Check mid-span section
            
            if sec_brittle_dict[i] == True:
                
                # Check demand moment against section ultimate moment stored in the capacity dict for such members
                
                if M_demand_mid > M_capacity_mid:
                    
                    rib_brittle_fail_dict = True # ultimate failure of system
                    
                    failed_rib_info.append((i, original_EI[i]))
                    
                    Flexure_fail_counter_mid +=1
                    
                    return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
                    
                   # print(f"Rib_{((i-num_points)/2)+1}: Compression-controlled failure at mid-span- system has failed")
                    
           
            # Check end span 1
            
            if sec_brittle_dict[i-num_points] == True:
                
                if M_demand_end > M_capacity_end_1:
                    
                    rib_brittle_fail_dict[i] = True
                    
                    failed_rib_info.append((i, original_EI[i]))
                    
                  #  print(f"Rib_{((i-num_points)/2)+1}: Compression-controlled failure at end 1- system has failed")
        
                    Flexure_fail_counter_end +=1
                    
                    return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
       
           
            # Check end span 2
        
            if sec_brittle_dict[i+num_points] == True:
            
               if M_demand_end > M_capacity_end_2:
                
                  rib_brittle_fail_dict[i] = True
                  
                  failed_rib_info.append((i, original_EI[i]))
                
                #  print(f"Rib_{((i-num_points)/2)+1}: Compression-controlled failure at end 2- system has failed")
                  
                  Flexure_fail_counter_end +=1
                  
                  return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
                  
           
            
            # Check end 1 if hinge already exists
            
            if plastic_hinge_end_1_dict[i] == True and rib_ductile_fail_dict[i] == False: 
                
              # print(f"Rib_{((i-num_points)/2)+1}: Plastic hinge already exists at end 1")
            
               M_redis_1 = M_demand_end - M_capacity_end_1 # moment in need of redistribution from each support
               
               if redis_limit * M_capacity_end_1 < M_redis_1:
                   
                  # print(f"Rib_{((i-num_points)/2)+1}: Moment redistribution limit exceeded for end 1- check if ultimate moment capacity is exceeded")
                   
                   M_redis_1 = redis_limit * M_capacity_end_1
                   
                   M_demand_end_1 = M_demand_end - M_redis_1
                   
                   if M_demand_end_1 > M_capacity_ult_end_1:
                       
                   #   print(f"Rib_{((i-num_points)/2)+1}: Ultimate moment capacity exceeded at end 1 - concrete crushing has occurred and system has failed")
                      
                      rib_brittle_fail_dict[i] = True
                      
                      failed_rib_info.append((i, original_EI[i]))
                      
                      Flexure_fail_counter_end +=1
                      
                      return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
                      
          
            
            # Check end 2 if hinge already exists
            
            if plastic_hinge_end_2_dict[i] == True and rib_ductile_fail_dict[i] == False: 
                
               # print(f"Rib_{((i-num_points)/2)+1}: Plastic hinge already exists at end 2")
                
                M_redis_2 = M_demand_end - M_capacity_end_2
                
                if redis_limit * M_capacity_end_2 < M_redis_2:
                    
                   # print(f"Rib_{((i-num_points)/2)+1}: Moment redistribution limit exceeded for end 2 - check if ultimate moment capacity is exceeded")
                    
                    M_demand_end_2 = M_demand_end - M_redis_2
                    
                    if M_demand_end_2 > M_capacity_ult_end_2:
                        
                     #  print(f"Rib_{((i-num_points)/2)+1}: Ultimate moment capacity exceeded at end 2 - concrete crushing has occurred and system has failed")
                       
                       rib_brittle_fail_dict[i] = True
                       
                       failed_rib_info.append((i, original_EI[i]))
                       
                       Flexure_fail_counter_end +=1
                       
                       return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
                       
                    
                       
            # Check if hinge exists at mid span
            
            if plastic_hinge_mid_dict[i] == True and rib_ductile_fail_dict[i] == False:
                
                #print(f"Rib_{((i-num_points)/2)+1}: Plastic hinge already exists at mid span")
                
                M_redis_mid = M_demand_mid - M_capacity_mid
                
                plastic_hinge_mid_dict[i] = True
                
                if redis_limit * M_capacity_mid < M_redis_mid:
                    
                  #  print(f"Rib_{((i-num_points)/2)+1}: Moment redistribution limit exceeded for mid-span - check if ultimate moment capacity is exceeded")
                    
                    M_demand_mid_ult = M_demand_mid - M_redis_mid
                    
                    if M_demand_mid_ult > M_capacity_ult_mid:
                        
                     #  print(f"Rib_{((i-num_points)/2)+1}: Ultimate moment capacity exceeded at mid-span - concrete crushing has occurred and system has failed")
                       
                       rib_brittle_fail_dict[i] = True
                       
                       failed_rib_info.append((i, original_EI[i]))
                       
                       Flexure_fail_counter_mid +=1
                       
                       return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
                       
                    
            # Work out new redistributed loads 
            
            # First case - plastic hinge at no ribs 
            
            if plastic_hinge_end_1_dict[i] == False and plastic_hinge_end_2_dict[i] == False and plastic_hinge_mid_dict[i] == False:
                
                M_demand_redis_end_1 = M_demand_end
                M_demand_redis_end_2 = M_demand_end
                M_demand_redis_mid = M_demand_mid
               
            # Second case - plastic hinge at end 1 only
            
            if plastic_hinge_end_1_dict[i] == True and plastic_hinge_end_2_dict[i] == False and plastic_hinge_mid_dict[i] == False:
                
                M_demand_redis_end_1 = 0 # already checked
                M_demand_redis_end_2 = M_demand_end 
                M_demand_redis_mid = M_demand_mid + M_redis_1 # all excess distributed to mid span
                
            # Third case - plastic hinge at end 2 only
            
            if plastic_hinge_end_1_dict[i] == False and plastic_hinge_end_2_dict[i] == True and plastic_hinge_mid_dict[i] == False:
                
                M_demand_redis_end_1 = M_demand_end 
                M_demand_redis_end_2 = 0 # already checked
                M_demand_redis_mid = M_demand_mid + M_redis_2 # all excess distributed to mid span
                
            # Fourth case - plastic hinge at mid-span only
            
            if plastic_hinge_end_1_dict[i] == False and plastic_hinge_end_2_dict[i] == False and plastic_hinge_mid_dict[i] == True:
                
                # Set up ratio based on moment capacity for sharing moment between two ends
                
                ratio_1 = M_capacity_end_1 / (M_capacity_end_1 + M_capacity_end_2) # proprtion of total excess moment from mid span redistributed to end 1
                ratio_2 = M_capacity_end_2 / (M_capacity_end_1 + M_capacity_end_2)
                
                M_demand_redis_end_1 = M_demand_end + M_redis_mid * ratio_1 
                M_demand_redis_end_2 = M_demand_end + M_redis_mid * ratio_2
                M_demand_redis_mid = 0 # already checked
             
            # Fifth case - plastic hinge at end 1 and end 2
            
            if plastic_hinge_end_1_dict[i] == True and plastic_hinge_end_2_dict[i] == True and plastic_hinge_mid_dict[i] == False: 
                
                M_demand_redis_end_1 = 0
                M_demand_redis_end_2 = 0
                M_demand_redis_mid = M_demand_mid + M_redis_1 + M_redis_2 # all moment has to go to mid span
             
            # Sixth case - plastic hinge at end 1 and mid-span
            
            if plastic_hinge_end_1_dict[i] == True and plastic_hinge_end_2_dict[i] == False and plastic_hinge_mid_dict[i] == True: 
                
                # Moment can only be redistributed practically from end to mid span - so only mid-span moments can be redistributed
                
                M_demand_redis_end_1 = M_demand_end + M_redis_1 # cannot be redistributed will lead to failure if ultimate capacity is surpassed
                M_demand_redis_end_2 = M_demand_end + M_redis_mid # takes moment redistributed from mid span only
                M_demand_redis_mid = 0
              
            # Seventh case - plastic hinge at end 2 and mid-span
            
            if plastic_hinge_end_1_dict[i] == False and plastic_hinge_end_2_dict[i] == True and plastic_hinge_mid_dict[i] == True: 
                
                M_demand_redis_end_1 = M_demand_end + M_redis_mid # only takes mid span moment
                M_demand_redis_end_2 = M_demand_end + M_redis_2 # cannot be practically redistributed
                M_demand_redis_mid = 0
              
            # Eigth case - plastic hinges at all sections
            
            if plastic_hinge_end_1_dict[i] == True and plastic_hinge_end_2_dict[i] == True and plastic_hinge_mid_dict[i] == True: 
               
                M_demand_redis_end_1 = 0
                M_demand_redis_end_2 = 0
                M_demand_redis_mid = 0
                
                rib_ductile_fail_dict[i] == True
                
            # Check for compression-controlled failure after redistribution
            
            # Check brittle sections for ultimate failure (they can have moment distributed to them but cannot redistribute moment themselves)
            
            # Check mid-span section
            
            if sec_brittle_dict[i] == True:
                
                # Check demand moment against section ultimate moment stored in the capacity dict for such members
                
                if M_demand_redis_mid > M_capacity_mid:
                    
                    rib_brittle_fail_dict = True # ultimate failure of system
                    
                    failed_rib_info.append((i, original_EI[i]))
                    
                    Flexure_fail_counter_mid +=1
                    
                   # print(f"Rib_{((i-num_points)/2)+1}: Compression-controlled failure at mid-span- system has failed")
                    
                    return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
                
            # Check end span 1
            
            if sec_brittle_dict[i-num_points] == True:
                
                if M_demand_redis_end_1 > M_capacity_end_1:
                    
                    rib_brittle_fail_dict[i] = True
                    
                    failed_rib_info.append((i, original_EI[i]))
                    
                    
                   # print(f"Rib_{((i-num_points)/2)+1}: Compression-controlled failure at end 1- system has failed")
                    
                    Flexure_fail_counter_end +=1
                    
                    return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
                  
        
            # Check end span 2
        
            if sec_brittle_dict[i+num_points] == True:
            
               if M_demand_redis_end_2 > M_capacity_end_2:
                
                  rib_brittle_fail_dict[i] = True
                  
                  failed_rib_info.append((i, original_EI[i]))
                
                 # print(f"Rib_{((i-num_points)/2)+1}: Compression-controlled failure at end 2- system has failed")
                  
                  Flexure_fail_counter_end +=1
                  
                  return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
                  
            
            
            # Check for new mid span hinge for next iteration
          
            if M_demand_redis_mid > M_capacity_mid and plastic_hinge_mid_dict[i] == False and rib_ductile_fail_dict[i] == False:
                
               # print(f"Rib_{((i-num_points)/2)+1}: New plastic hinge has formed at mid-span")
                
                plastic_hinge_mid_dict[i] = True
                
                if M_demand_redis_mid > M_capacity_ult_mid:
                    
                  #  print(f"Rib_{((i-num_points)/2)+1}: Moment redistribution to mid-span has caused ultimate failure")
                    
                    rib_brittle_fail_dict[i] = True
                    
                    failed_rib_info.append((i, original_EI[i]))
                    
                    Flexure_fail_counter_mid +=1
                    
                    return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
                    
               
                
            # Check for new end 1 hinge
            
            if M_demand_redis_end_1 > M_capacity_end_1 and plastic_hinge_end_1_dict[i] == False and rib_ductile_fail_dict[i] == False:
                
               # print(f"Rib_{((i-num_points)/2)+1}: New plastic hinge has formed at end 1")
                
                plastic_hinge_end_1_dict[i] = True
                
                if M_demand_redis_end_1 > M_capacity_ult_end_1:
                    
                   # print(f"Rib_{((i-num_points)/2)+1}: Moment redistribution to end 1 has caused ultimate failure")
                    
                    rib_brittle_fail_dict[i] = True
                    
                    failed_rib_info.append((i, original_EI[i]))
                    
                    Flexure_fail_counter_end +=1
                    
                    return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
           
             
            # Check for new end 2 hinge
            
            if M_demand_redis_end_2 > M_capacity_end_2 and plastic_hinge_end_2_dict[i] == False and rib_ductile_fail_dict[i] == False:
                
               # print(f"Rib_{((i-num_points)/2)+1}: New plastic hinge has formed at end 2")
                
                plastic_hinge_end_2_dict[i] = True
                
                if M_demand_redis_end_2 > M_capacity_ult_end_2:
                    
                  #  print(f"Rib_{((i-num_points)/2)+1}: Moment redistribution to end 2 has caused ultimate failure")
                    
                    rib_brittle_fail_dict[i] = True
                    
                    failed_rib_info.append((i, original_EI[i]))
                    
                    Flexure_fail_counter_end +=1
                    
                    return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
                    
                    
                
            # Adjust stiffness for next load application - if number of hinges has increased - these hinges can be evaluated next time around
            
            # First case - plastic hinge at no ribs 

            # no stiffness adjustement required
                
            # Second case - plastic hinge at end 1 only
            
            if plastic_hinge_end_1_dict[i] == True and plastic_hinge_end_2_dict[i] == False and plastic_hinge_mid_dict[i] == False:
     
                # Consider propped cantilever case
                
                loc = 0.5
             
                k = 48 * EI / (((span ** 3) - 3 * span * ((loc * span) ** 2) + 2 * ((loc * span) ** 3)))
                
                k_dict[i] = k
                
              #  print(f"Rib_{((i-num_points)/2)+1}: Hinge at end 1 - update to propped cantilever stiffness")
                
            # Third case - plastic hinge at end 2 only
            
            if plastic_hinge_end_1_dict[i] == False and plastic_hinge_end_2_dict[i] == True and plastic_hinge_mid_dict[i] == False:
     
                loc = 0.5
             
                k = 48 * EI / (((span ** 3) - 3 * span * ((loc * span) ** 2) + 2 * ((loc * span) ** 3)))
                
                k_dict[i] = k
                
               # print(f"Rib_{((i-num_points)/2)+1}: Hinge at end 2 - update to propped cantilever stiffness")
                
            # Fourth case - plastic hinge at mid-span only
            
            if plastic_hinge_end_1_dict[i] == False and plastic_hinge_end_2_dict[i] == False and plastic_hinge_mid_dict[i] == True:

                k = 64 * EI / (span ** 3) # cantilever end span deflection but considered for half span
           
                k_dict[i] = k
                
              #  print(f"Rib_{((i-num_points)/2)+1}: Hinge at mid-span - update to end span cantilever stiffness")
                
            # Fifth case - plastic hinge at end 1 and end 2
            
            if plastic_hinge_end_1_dict[i] == True and plastic_hinge_end_2_dict[i] == True and plastic_hinge_mid_dict[i] == False: 
          
                k = (384 / 5) * EI / (span ** 3)
             
                k_dict[i] = k
            
               # print(f"Rib_{((i-num_points)/2)+1}: Two plastic hinges - update to simply supported stiffness")
                
            # Sixth case - plastic hinge at end 1 and mid-span
            
            if plastic_hinge_end_1_dict[i] == True and plastic_hinge_end_2_dict[i] == False and plastic_hinge_mid_dict[i] == True:
           
                k = (384 / 5) * EI / (span ** 3)
              
                k_dict[i] = k
             
                print(f"Rib_{((i-num_points)/2)+1}: Two plastic hinges - update to simply supported stiffness")
                
            # Seventh case - plastic hinge at end 2 and mid-span
            
            if plastic_hinge_end_1_dict[i] == False and plastic_hinge_end_2_dict[i] == True and plastic_hinge_mid_dict[i] == True: 
      
                k = (384 / 5) * EI / (span ** 3)
             
                k_dict[i] = k
             
                #print(f"Rib_{((i-num_points)/2)+1}: Two plastic hinges - update to simply supported stiffness")
                
            # Eigth case - plastic hinges at all sections
            
            if plastic_hinge_end_1_dict[i] == True and plastic_hinge_end_2_dict[i] == True and plastic_hinge_mid_dict[i] == True: 
                
               # print(f"Rib_{((i-num_points)/2)+1}: Three plastic hinges formed - reduce to residual stiffness and rely on adjacent ribs to avoid collapse")
                
                rib_ductile_fail_dict == True
                
                k = residual_factor * 384 * EI / (span ** 3)
                
                k_dict[i] = k
                
                collapse_counter += 1
              
                
    deflection_dict[i] = max(deflection_max_dict)
    
    return System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
                
             # Now reapply same load to evaluate these hinges
                    
def render_model(model_1):
    
    return
 
    renderer = Renderer(model_1)
    renderer.annotation_size = 50
    renderer.deformed_shape = True
    renderer.deformed_scale = 100
    renderer.render_loads = False
    renderer.combo_name = 'combo'
    renderer.labels = False
    renderer.screenshot("model_1.png",interact=True)
    renderer.render_model()
                     
def iterate_full_load(System_ULS_fail_counter, k_dict, rib_ductile_fail_dict, rib_brittle_fail_dict, rib_status_dict, E_c_dict, I_dict, scaled_node_load_dict, strain_c_ult_dict, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter):
    
    # Apply full load
    
    #print('Full load applied, iterate to find stable state')
    
    scaled_node_load_dict = [node_load for node_load in node_load_dict]
    
    model_1 = FEA(E_c_dict, I_dict, k_dict, scaled_node_load_dict)
    
    model_1.analyze()
    
    # Process output
  
    System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter = rib_output(k_dict, rib_ductile_fail_dict, rib_brittle_fail_dict, rib_status_dict, model_1, strain_c_ult_dict, System_ULS_fail_counter, x_dict, M_capacity_dict, EI_dict, d_dict, strain_s_ult_dict, plastic_hinge_end_1_dict, plastic_hinge_end_2_dict, plastic_hinge_mid_dict, v_capacity_dict, b_dict, z_dict, h_f_dict, tau_capacity_dict, f_c_dict, M_capacity_ult_dict, deflection_max_dict, sec_brittle_dict, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter)
                
    # Check for brittle failure
    
    if sum(rib_brittle_fail_dict) > 0:
      
        System_ULS_fail_counter += 1
       
      #  print(f'{sum(rib_brittle_fail_dict)} ultimate failures have occurred, system has failed')
        
        return System_ULS_fail_counter,Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter

        
    else: # if no brittle failures occur
    
         num_failures = sum(rib_ductile_fail_dict)
         
        # print(f'{num_failures} plastic hinges at mid span')
    
         # Re run with adjusted stiffness but same load

         while True: # number of failures still increasing, loop infinfitely until break out where number of failures stop increasing
         
               num_failures_org = num_failures
         
               # Re run with adjusted stiffness but same load
                   
               # Input to FEA
                   
               model_1 = FEA(E_c_dict, I_dict, k_dict, scaled_node_load_dict)
                   
               # Conduct FEA
                       
               model_1.analyze()
               
               # Render the model
               
               render_model(model_1)
         
               # Process output
               
               System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter = rib_output(k_dict, rib_ductile_fail_dict, rib_brittle_fail_dict, rib_status_dict, model_1, strain_c_ult_dict, System_ULS_fail_counter, x_dict, M_capacity_dict, EI_dict, d_dict, strain_s_ult_dict, plastic_hinge_end_1_dict, plastic_hinge_end_2_dict, plastic_hinge_mid_dict, v_capacity_dict, b_dict, z_dict, h_f_dict, tau_capacity_dict, f_c_dict, M_capacity_ult_dict, deflection_max_dict, sec_brittle_dict, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter)
        
               # Check for brittle and hence system failure of any ribs
               
               if sum(rib_brittle_fail_dict) > 0:
                   
                   System_ULS_fail_counter += 1
                   
                  # print(f'{sum(rib_brittle_fail_dict)} ultimate failures have occurred, system has failed')
                   
                   return System_ULS_fail_counter,Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
               
               # if no brittle failures occur
                   
               # First check number of total failures
                   
               new_failures = sum(rib_ductile_fail_dict) 
                   
              # print(f'{new_failures} plastic hinges have formed at mid span')
               
               if new_failures == num_failures_org:
                     
                  break
              
               # Update for next iteration if the new failure count has changed
               
               num_failures = new_failures
               
    return System_ULS_fail_counter,Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
               
def iterate_scale_factor(iter_counter, System_ULS_fail_counter, k_dict, rib_ductile_fail_dict, rib_brittle_fail_dict, rib_status_dict, E_c_dict, I_dict, scaled_node_load_dict, strain_c_ult_dict, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter):
    
    # Scale factor is not one so it must be less than one
    
    # Calculate load increment of ten parts to full load
    
    load_inc_counter = 0
    
    # Establish current scale factor
    
   # print(f'Initial scale factor is {scale_factor}')

    scaled_node_load_dict = [scale_factor * node_load for node_load in node_load_dict]
    
    model_1 = FEA(E_c_dict, I_dict, k_dict, scaled_node_load_dict)
    
    model_1.analyze()
    
    # Process output

    System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter = rib_output(k_dict, rib_ductile_fail_dict, rib_brittle_fail_dict, rib_status_dict, model_1, strain_c_ult_dict, System_ULS_fail_counter, x_dict, M_capacity_dict, EI_dict, d_dict, strain_s_ult_dict, plastic_hinge_end_1_dict, plastic_hinge_end_2_dict, plastic_hinge_mid_dict, v_capacity_dict, b_dict, z_dict, h_f_dict, tau_capacity_dict, f_c_dict, M_capacity_ult_dict, deflection_max_dict, sec_brittle_dict, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter)
    
    # Check for brittle and hence system failure of any ribs
    
    if sum(rib_brittle_fail_dict) > 0:
        
        System_ULS_fail_counter += 1
        
       # print(f'{sum(rib_brittle_fail_dict)} ultimate failures have occurred, system has failed')
    
    else:  # If no brittle failures occur
    
        # First check number of total failures
        
        num_failures = sum(rib_ductile_fail_dict) 
        
      #  print(f'{num_failures} plastic hinges have formed at mid span')
        
        render_model(model_1)

        # Iterate
        
        while load_inc_counter < num_inc:  
        
            num_failures_org = num_failures  # Track the last failure count
        
            # Re-run FEA with adjusted stiffness but same load
            
            model_1 = FEA(E_c_dict, I_dict, k_dict, scaled_node_load_dict)
            
            model_1.analyze()
            
            # Render the model
            
            render_model(model_1)

            # Process output
           
            System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter = rib_output(k_dict, rib_ductile_fail_dict, rib_brittle_fail_dict, rib_status_dict, model_1, strain_c_ult_dict, System_ULS_fail_counter, x_dict, M_capacity_dict, EI_dict, d_dict, strain_s_ult_dict, plastic_hinge_end_1_dict, plastic_hinge_end_2_dict, plastic_hinge_mid_dict, v_capacity_dict, b_dict, z_dict, h_f_dict, tau_capacity_dict, f_c_dict, M_capacity_ult_dict, deflection_max_dict, sec_brittle_dict, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter)
            
            # Check for brittle failure
            
            if sum(rib_brittle_fail_dict) > 0:
              
                System_ULS_fail_counter += 1
               
               # print(f'{sum(rib_brittle_fail_dict)} ultimate failures have occurred, system has failed')
              
                return System_ULS_fail_counter,Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter
              # Exit immediately if brittle failure occurs
        
            # Check total failures
           
            num_failures = sum(rib_ductile_fail_dict) 
           
           # print(f'{num_failures} plastic hinges have formed at mid span')
        
            # If no additional failures occur, increase load
        
            if num_failures == num_failures_org:
               
                print(f'Steady state reached with {num_failures} plastic hinges at mid span, increase load')
                      
                # Update load dictionary
                
                load_inc_counter += 1
               
                scaled_node_load_dict = [
                   scale_factor * node_load + 
                    load_inc_counter * ((node_load - scale_factor * node_load) / num_inc) 
                    for node_load in node_load_dict
                ]
        
            if load_inc_counter == num_inc:
                
              #  print("Max load increments reached. Beginning full load iterations.")
                
                System_ULS_fail_counter,Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter = iterate_full_load(System_ULS_fail_counter, k_dict, rib_ductile_fail_dict, rib_brittle_fail_dict, rib_status_dict, E_c_dict, I_dict, scaled_node_load_dict, strain_c_ult_dict, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter)
            
    return System_ULS_fail_counter,Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter

# Simulation parameters

num_sims = 100000
seed = 1
np.random.seed(seed)
scale_factor = 1.0 # initial scale factor defined for load increment (all at once)
num_inc = 1
reduction_factor = 1

# Model parameters

d_load = 0.0015 # typical superimposed dead load
span = 6000
pot = 300
num_ribs = 10
b_nom = 125 # nominal rib width from Southampton Way roof extension report by MBP
h_nom = 250 # nominal rib depth from Southampton Way roof extension report by MBP
h_f_nom = 75 # nominal topping slab depth from Southampton Way roof extension report by MBP
c = 25 # nominal concrete cover according to CP 114
c_top_dev_mu = 10 # mean deviation from design top concrete cover, highly variable in a range of 5 to 16 according to JCSS PMC 3.2.2 Table 1 for UK so should be defined appropriately by user
c_bot_dev_mu = -7 # mean deviation from design bottom concrete cover, highly variable in a range of 7 to 23 according to JCSS PMC 3.2.2 Table 1 for UK so should be defined appropriately by user
A_s_nom = 628.3 # assume two 20mm bars for now but this can change with discrete random variable for steel not included here
rib_spacing = 425 # distance c/c -assumed not to be influenced by geometric random variables as it is large
num_slabs = num_ribs - 1
num_points = num_ribs + num_slabs
num_grid = num_points * 3
point_spacing = rib_spacing / 2 # distance
l_0_factor = 0.5 # isolated encastre case
residual_factor = 0.05 # True hinge
redis_limit = 0.30 # 30% limit on moment redistribution within each rib
conf_factor = 0 # derived from relation between longitudinal stress and confining stress
cal_factor = 1

# Concrete deterministic parameters

poissons_ratio = 0.2 # standard assumption for concrete in elastic range
t_load = 28 # age at time of loading in days assumed for now
tau_yr = 70 # duration of loading in years also assumed
tau = tau_yr * 365.25
RH = 50 # percentage for typical environment described in Eurocode 2
D_lower = 20 # smallest aggregate size in top sieve, assume 20 for now has a notable impact on concrete shear capacity

# Creep factor definition

# Assuming CEM 1 ordinary portland cement, CN

sc_creep = 0.5 # low strength class less than 35 MPa
alpha_sc = 0 # 0 for normal concrete

# Consider basic creep coefficient and drying creep coefficient

creep_factor_cov = 0.3 # defined in eurocode 2

# Reinforcement

# Bottom reinforcement in the ribs

# CP114 specifies a minim bar diameter of 16mm in the ribs

# BS785:1967 Hot Rolled Steel Bars

# Assume use of mild steel

# Yield stress for bars between 9.5 and 38.1 mm  diameter (coverted fron inches)

S_nom = 248.11 # converted to MPa from kgf/m^2

# Assumed placement of bottom steel in ribs at mid span

# Two 16 mm bars or one 20 mm bar, assume 90% and 10% chance, respectively

num_bars_mid_p1 = 2
num_bars_mid_p2 = 1

bar_diameter_mid_p1 = 16
bar_diameter_mid_p2 = 20

mid_p1 = 0.8
mid_p2 = 0.2

# Assumed placement of steel fabric in topping slab BS1221:1945 or mild round steel bars

# Go with roud mild steel bars

# 8mm bars at 200mm c/c, 6 in one metre section

bar_diameter_slab = 8
spacing_slab = 300
num_bars_slab = int((span / spacing_slab) + 1)

# Assumed end bars

num_bars_end = 2

bar_diameter_end_p1 = 20
bar_diameter_end_p2 = 16

end_p1 = 0.8
end_p2 = 0.2
            
# Concrete material property definition

# Determinstic parameter alpha

# Calculate alpha_1 which is influenced by the duration of loading

# Define input parameters

a_tau = 0.04 # given in JCSS PMC Section 3.1.1

alpha_3_inf = 0.8 # given as approximately 0.8 in JCSS PMC Section 3.1.1

alpha_1 = alpha_3_inf + (1 - alpha_3_inf) * np.exp(-a_tau * tau) # Eq. 5a will give 0.8 for almost all cases so could be taken deterministicially

# Caclulate alpha_2 which is influenced by the age at time of loading

# Define input parameters for normal conditions as given in JCSS PMC Section 3.1.1

a = 0.6
b = 0.12

alpha_2 = a + b * np.log(t_load) # Eq. 5b

# calculate alpha

alpha = alpha_1 * alpha_2

# Define lambda accounting for systematic variation between in situ compressive strength and standard compression test results

# log normal random variable but JCSS allows it to be take determinstically 

lambda_fc = 0.96

# Define random variable for concrete strength

# Set up a random field for conrete strength variability within one member, represented by standard normal variable U

# Correlation parameters

rho = 0.5 # correlation factor given in JCSS PMC for Eq.17
d_c = 6000 # converted to mm, correlation length given in JCSS PMC for Eq.17

# Generate empty correllation matrix for the point grid

correlation_matrix = np.zeros((num_grid, num_grid))

# Define spacing of points

spacing_y = span / 2
spacing_x = point_spacing

# Define number of points

num_y = 3
num_x = num_points

# Create coordinate sets

x_coords = np.arange(num_x) * spacing_x
y_coords = np.arange(num_y) * spacing_y

# Generate grid

xx, yy = np.meshgrid(x_coords, y_coords)
node_matrix = np.column_stack((xx.ravel(), yy.ravel()))

# Populate the correlation matrix using Eq. 17 from JCSS PMC

for i in range(num_grid):
 
    for j in range(num_grid):
        
        r_ij = np.sqrt((node_matrix[i, 0] - node_matrix[j, 0])**2 + (node_matrix[i, 1] - node_matrix[j, 1])**2)

        correlation_matrix[i, j] = rho + (1 - rho) * np.exp((-(r_ij) ** 2)/ ((d_c) ** 2))

# Carry out eigen decompostition of the correlation matrix to find correlation_matrix = L * L.T
 
eigen_values, eigen_vectors = np.linalg.eigh(correlation_matrix)
 
# Add a tolerance threshold to account for computer inaccuracy
 
tolerance = 1e-15

# Set any eigen values below tolerance to 0

eigen_values[eigen_values < tolerance] = 0

# Calculate L where L = T @ D_sqrt
 
# Assign matrix of eigen vectors to variable T
 
T = eigen_vectors
 
# Generate the root of the diagonal matrix of eigen values
 
D_sqrt = np.diag(np.sqrt(eigen_values))
 
# Caclulate the matrix L
 
L = T @ D_sqrt # L can then be applied ot the matrix of uncorrelated standard normals for each simulation

# Define Student's T distribution for basic concrete compression strength f_co from JCSS PMC 3.1.5

# Assume concrete strength class C20/25 as its properties best corresponds with 1:2:4 ready mix in CP114

# Prior parameters if no specific information is available from JCSS PMC Table 3.1.2

# As C20/25 not featured use updated parameters from Fieri et al, 2023

# Parameter definitions from in Fieri et al, 2023

m = 3.525 # logarithmic mean 
n = 3 # equivalent sample size
s = 0.13 # empircial standard deviation of an equivalent sample of size v + 1
v = 10 # degrees of freedom of the Student's T distribution

# Define log normal random variable Y_1j accounting for special placing, curing and hardening conditions as defined in JCSS PMC for compression strength

# Normal space parameters as per JCSS PMC Table 3.1.1

Y_1j_mean = 1
Y_1j_cov = 0.06
Y_1j_std = Y_1j_mean * Y_1j_cov

# Convert to log space parameters using widely known equations 

Y_1j_ln_var = np.log(1.0 + (Y_1j_std ** 2) / (Y_1j_mean ** 2))
Y_1j_ln_std = np.sqrt(Y_1j_ln_var)
Y_1j_ln_mu = np.log(Y_1j_mean) - 0.5 * Y_1j_ln_var

# Define Y_2j for tensile strength

# Normal space parameters as per JCSS PMC Table 3.1.1

Y_2j_mean = 1
Y_2j_cov = 0.3
Y_2j_std = Y_2j_mean * Y_2j_cov

# Convert to log space parameters using widely known equations 

Y_2j_ln_var = np.log(1.0 + (Y_2j_std ** 2) / (Y_2j_mean ** 2))
Y_2j_ln_std = np.sqrt(Y_2j_ln_var)
Y_2j_ln_mu = np.log(Y_2j_mean) - 0.5 * Y_2j_ln_var

# Define Y_3j for the Young's Modulus

# Normal space parameters as per JCSS PMC Table 3.1.1

Y_3j_mean = 1
Y_3j_cov = 0.15
Y_3j_std = Y_3j_mean * Y_3j_cov

# Convert to log space parameters using widely known equations 

Y_3j_ln_var = np.log(1.0 + (Y_3j_std ** 2) / (Y_3j_mean ** 2))
Y_3j_ln_std = np.sqrt(Y_3j_ln_var)
Y_3j_ln_mu = np.log(Y_3j_mean) - 0.5 * Y_3j_ln_var

# Define y_4j for ultimate strain

# Normal space parameters as per JCSS PMC Table 3.1.1

Y_4j_mean = 1
Y_4j_cov = 0.15
Y_4j_std = Y_4j_mean * Y_4j_cov

# Convert to log space parameters using widely known equations 

Y_4j_ln_var = np.log(1.0 + (Y_4j_std ** 2) / (Y_4j_mean ** 2))
Y_4j_ln_std = np.sqrt(Y_4j_ln_var)
Y_4j_ln_mu = np.log(Y_4j_mean) - 0.5 * Y_4j_ln_var

# Steel material property defintion, only varies between each bar, not longitudinally along bars

# Define determinsitic values

poissons_ratio_s = 0.3 # widely accepted assumption

# Define normally distributed random variable for older steel yield strength fy according to Wisniewski et al, 2009 for older steels

lambda_s = 1.21 # bias of mean relative to nominal strength for older steels as given in Wisniewski et al, 2009

f_y_mean = S_nom * lambda_s # occurs at 2% strain

f_y_cov = 0.1 # for older steel 10%

f_y_std = f_y_mean * f_y_cov

# Define normally distributed random variable for steel yield strength fy according to Wisniewski et al, 2009 for older steels

# 15-20% higher than yield strength, take mid-point

f_ult_mean = f_y_mean * 1.175 # 17.5 % increase

f_ult_cov = 0.1 # again taken as 10% for older steels

f_ult_std = f_ult_mean * f_ult_cov

# Define normally distributed random variable for ultimate strain

strain_s_ult_mean = 0.1 # 10%

strain_s_ult_cov = 0.15 # 15%

strain_s_ult_std = strain_s_ult_mean * strain_s_ult_cov

# Define normally distributed random variable for Young's Modulus as per Wisniewski et al, 2009

# Shows no correlation with other properties based on Wisniewski et al, 2009 Figure 7

E_s_mean = 202000 # converted from GPa to MPa

E_s_cov = 0.04 # 4%

E_s_std = E_s_mean * E_s_cov

# Define normally distributed steel area random variable as per JCSS PMC Table 3.2.2

A_s_mean = A_s_nom # mean is equal to nominal steel area 

A_s_cov = 0.02 # 2%

A_s_std = A_s_mean * A_s_cov

# Set up correlation framework between material properties as outlined in JCSS PMC Table 3.2.2

# Consider matrix of f_y, f_ult, strain_s_ult, A_s

# Create empty correlation matrix

steel_correlation_matrix = np.zeros((4, 4)) # Four material/sectional properties are considered

# Define correlation coefficient between each variable from JCSS PMC Table 3.2.2

rho_A_s_f_y = 0.5
rho_A_s_f_ult = 0.35
rho_A_s_strain_s_ult = 0
rho_f_y_f_ult = 0.85
rho_f_y_strain_s_ult = - 0.5
rho_f_ult_strain_s_ult = - 0.55

# Populate correlation matrix

# Define perfect positive correlation between a property and itself

steel_correlation_matrix[0,0] = 1 # avoid indexing error
steel_correlation_matrix[1,1] = 1
steel_correlation_matrix[2,2] = 1
steel_correlation_matrix[3,3] = 1

# Add correlation coefficients

steel_correlation_matrix[0,1] = rho_A_s_f_y
steel_correlation_matrix[0,2] = rho_A_s_f_ult
steel_correlation_matrix[1,0] = rho_A_s_f_y
steel_correlation_matrix[2,0] = rho_A_s_f_ult
steel_correlation_matrix[1,2] = rho_f_y_f_ult 
steel_correlation_matrix[1,3] = rho_f_y_strain_s_ult
steel_correlation_matrix[2,1] = rho_f_y_f_ult 
steel_correlation_matrix[3,1] = rho_f_y_strain_s_ult
steel_correlation_matrix[2,3] = rho_f_ult_strain_s_ult
steel_correlation_matrix[3,2] = rho_f_ult_strain_s_ult

# Conduct eigen decomposition as done for the correlation matrix of concrete strength

eigen_values_s, eigen_vectors_s = np.linalg.eigh(steel_correlation_matrix)

# Set any eigen values below tolerance to 0

eigen_values_s[eigen_values_s < tolerance] = 0

# Calculate L where L = T @ D_sqrt
 
# Assign matrix of eigen vectors to variable T
 
T_s = eigen_vectors_s
 
# Generate the root of the diagonal matrix of eigen values
 
D_sqrt_s = np.diag(np.sqrt(eigen_values_s))
 
# Caclulate the matrix L for steel material poperties
 
L_s = T_s @ D_sqrt_s # L can then be applied to the vector of properties created during each simulation for each point

# Dimensions

# According to JCSS no strong correlation exists between horizontal and vertical dimensions so can be ignored

# External dimensions of concrete components

# Nominal dimensions defined in model parameters section

# Define normnally distributed random variables for deviation from nominal dimensions as per JCSS PMC Section 3.10.1

h_dev_mu = 0.003 * h_nom # but no less than 0 abnd no more than 3 mm

if h_dev_mu < 0:
    
    h_dev_mu = 0

if h_dev_mu > 3:
    
    h_dev_mu = 3
    
b_dev_mu = 0.003 * b_nom # but no less than 0 abnd no more than 3 mm

if b_dev_mu < 0:
    
    b_dev_mu = 0

if b_dev_mu > 3:
    
    b_dev_mu = 3
    
h_f_dev_mu = 0.003 * h_f_nom # but no less than 0 abnd no more than 3 mm

if h_f_dev_mu < 0:
    
    h_f_dev_mu = 0

if h_f_dev_mu > 3:
    
    h_f_dev_mu = 3

h_dev_std = 4 + 0.006 * h_nom # but no more than 10 mm

if h_dev_std > 10:
    
    h_dev_std = 10
    
b_dev_std = 4 + 0.006 * b_nom # but no more than 10 mm

if b_dev_std > 10:
    
    b_dev_std = 10
    
h_f_dev_std = 4 + 0.006 * h_f_nom # but no more than 10 mm

if h_f_dev_std > 10:
    
    h_f_dev_std = 10
    
# Concrete cover for top steel

c_nom = 25

c_top_dev_std = 11 # JCSS PMC 3.3.3 Table 1 gives a range of 6 to 16 take mid point of 11

# Concrete cover for bottom steel

# Mean is highly variable so defined in model parameters

c_bot_dev_std = 3.5 # given in JCSS PMC 3.2.2 Table 1 as 3-4 take mid point of 3.5

# Loading

# Sustained live load

# Parameters from JCSS PMC Table 2.2.1

# Load intensity parameters W

m_q = 0.5 # for residential sustained live load
sigma_V = 0.3 # standard deviation of variable V
sigma_U_q = 0.6

var_V = sigma_V ** 2
var_U_q = sigma_U_q ** 2

# Mean and variance of gamma distribution

E_q = m_q
var_q = var_V + var_U_q
mew_U_q = E_q / var_q
k_q = (E_q ** 2) / var_q

# Parameters for possion process

T_yr = 1 # extended service life of floor

lambda_q = 1 / 5

# Intermittent live load

m_p = 0.2
sigma_U_p = 0.4
var_U_p = sigma_U_p ** 2
var_p = var_V + var_U_p
E_p = m_p
mew_U_p = E_p / var_p
k_p = (E_p ** 2) / var_p

# Parameters for poisson process

v_p = 1 / 0.3


# Set up a random field for conrete strength variability within one member, represented by standard normal variable U

correlation_matrix_q = np.zeros((num_grid, num_grid))

# Populate the correlation matrix for W_q using Eq. 17 from JCSS PMC

for i in range(num_grid):
 
    for j in range(num_grid):
        
        r_ij = np.sqrt((node_matrix[i, 0] - node_matrix[j, 0])**2 + (node_matrix[i, 1] - node_matrix[j, 1])**2)

        correlation_matrix_q[i, j] = (np.exp((-(r_ij) ** 2)/ ((d_c) ** 2)) * var_U_q + var_V) / var_q 
        
correlation_matrix_p = np.zeros((num_grid, num_grid))

# Populate the correlation matrix for W_q using Eq. 17 from JCSS PMC

for i in range(num_grid):
 
    for j in range(num_grid):
        
        r_ij = np.sqrt((node_matrix[i, 0] - node_matrix[j, 0])**2 + (node_matrix[i, 1] - node_matrix[j, 1])**2)

        correlation_matrix_p[i, j] = (np.exp((-(r_ij) ** 2)/ ((d_c) ** 2)) * var_U_p + var_V) / var_p 

# Carry out eigen decompostition of the correlation matrix to find correlation_matrix = L * L.T
 
# Add a tolerance threshold to account for computer inaccuracy

tolerance = 1e-15

# Carry out eigen decomposition of the correlation matrix

eigen_values_q, eigen_vectors_q = np.linalg.eigh(correlation_matrix_q)

# Set any eigen values below tolerance to 0

eigen_values_q[eigen_values_q < tolerance] = 0

# Calculate L where L = T @ D_sqrt

T_q = eigen_vectors_q
D_q_sqrt = np.diag(np.sqrt(eigen_values_q))
L_q = T_q @ D_q_sqrt

# Carry out eigen decomposition of the correlation matrix

eigen_values_p, eigen_vectors_p = np.linalg.eigh(correlation_matrix_p)

# Set any eigen values below tolerance to 0

eigen_values_p[eigen_values_p < tolerance] = 0

# Calculate L where L = T @ D_sqrt

T_p = eigen_vectors_p
D_p_sqrt = np.diag(np.sqrt(eigen_values_p))
L_p = T_p @ D_p_sqrt


# Run Monte Carlo Simulation

System_SLS_fail_counter = 0
System_ULS_fail_counter = 0
Flexure_fail_counter_mid = 0
Flexure_fail_counter_end = 0
Vshear_counter = 0
Lshear_counter = 0
collapse_counter = 0

failed_rib_info = []

# Start timing
import time
start_time = time.time()


for sim in range(num_sims):
    
    rib_ductile_fail_dict = [False] * num_grid
    rib_brittle_fail_dict = [False] * num_grid
    rib_ductile_fail_dict_org = rib_ductile_fail_dict.copy()
    rib_brittle_fail_dict_org = rib_brittle_fail_dict.copy()
    sec_brittle_dict = [False] * num_grid

    # Set up dictionaries for key parameters and capacities
    
    # Capacity stores
    
    M_capacity_dict = [0] * num_grid
    M_capacity_ult_dict = [0] * num_grid
    v_capacity_dict = [0] * num_grid
    tau_capacity_dict = [0] * num_grid 
    node_load_max_dict = [0] * num_grid
    deflection_max_dict = [0] * num_grid
    original_EI = [0]*num_grid
    
    # Load stores
    
    ext_node_load_dict = [0] * num_grid # live and dead load combined
    sw_dict = [0] * num_grid
    node_load_dict = [0] * num_grid
    sw_node_load_dict = [0] * num_grid
    
    # Concrete material properties
    
    f_c_dict = [0] * num_grid
    E_c_dict = [0] * num_grid
    strain_c_ult_dict = [0] * num_grid
    creep_factor_dict = [0] * num_grid
    Y_3j_dict = [0] * num_grid
  
    # Load dict
    
    q_dict = [0] * num_grid
    p_dict = [0] * num_grid
 
    # Steel material properties and cover

    f_y_dict = []
    E_s_dict = []
    A_s_dict = []
    A_s_total_dict = [0] * num_grid
    c_bot_dict = [0] * num_grid
    c_top_dict = [0] * num_grid
    strain_s_ult_dict = []
    strain_s_y_dict = []
   
    # Dimensions of concrete parts 
    
    h_dict = [0] * num_grid
    b_dict = [0] * num_grid
    h_f_dict = [0] * num_grid
    b_eff_dict = [0] * num_grid

    # Concrete section properties
    
    z_dict = [0] * num_grid
    x_dict = [0] * num_grid
    I_dict = [0] * num_grid
    d_dict = [0] * num_grid
    rib_status_dict = [0] * num_grid
    EI_dict = [0] * num_grid
    EI_eq_dict = [0] * num_grid

    # Stiffness properties at mid span
    
    k_dict = [0] * num_grid # carries updated stiffness values
    k_dict_org = [0] * num_grid # keeps original stiffnesses until a successful redefintion takes place for a subsequent step
    
    # Output stores
    
    M_demand_end_dict = [0] * num_grid
    M_demand_mid_dict = [0] * num_grid
    V_demand_dict = [0] * num_grid
    deflection_dict = [0] * num_grid
    plastic_hinge_end_1_dict = [False] * num_grid
    plastic_hinge_end_2_dict = [False] * num_grid
    plastic_hinge_mid_dict = [False] * num_grid
    
    current_seed = seed + sim # record seed number for printing
    
    print(f"\n--- Simulation {sim+1}, Seed: {current_seed} ---\n")
    
    # Generate random variables
    
    # Generate vector of uncorrelated standard normals for U
    
    U_uncorr_vector = np.random.randn(num_grid)
    
    # Apply matrix L to the uncorrelated vector to generate the correlated vector of standard normals U
    
    U_vector = L @ U_uncorr_vector
    
    # Generate Y_1j for concrete compressive strength
    
    Y_1j = np.random.lognormal(Y_1j_ln_mu, Y_1j_ln_std)
    
    # Generate Y_2j for concrete tensile strength
    
    Y_2j = np.random.lognormal(Y_2j_ln_mu, Y_2j_ln_std)
    
    # Generate Y_3j for concrete Young's Modulus
    
    Y_3j = np.random.lognormal(Y_3j_ln_mu, Y_3j_ln_std)
    
    # Generate Y_4j for concrete ultimate strain
    
    Y_4j = np.random.lognormal(Y_4j_ln_mu, Y_4j_ln_std)
    
    # Calculate standard deviation according to Student's T distibution
    
    # Take a single sample from Student's T distribution
    
    t_v = t.rvs(v)
    
    sigma_j = t_v * s * (((1 + (1 / n)) ** 0.5))
    
    # Generate sustained load events
    
    # Consider number of sustained states with the exponential distribution
    
    T_total = T_yr      # 10-year reference period

  
    T_current = 0
    
    durations_q = []
    starts_q = []
    ends_q = []

    while T_current < T_total:
        
        duration_q = np.random.exponential(scale=1/lambda_q)
        remaining_time = T_total - T_current
        actual_duration = min(duration_q, remaining_time)
        
        durations_q.append(actual_duration)
        starts_q.append(T_current)
        T_current += actual_duration
        ends_q.append(T_current)
    
    N_q_events = len(durations_q)

    W_q_event_matrix = np.zeros((N_q_events, num_grid))

    for i in range(N_q_events):
        # 1. Generate independent standard normal vector
       
        W_q_uncorr_vector_i = np.random.normal(0, 1, size=num_grid)

        # 2. Apply transformation to get correlated standard normals
       
        W_q_vector_norm_i = L_q @ W_q_uncorr_vector_i

        # 3. Transform to gamma-distributed sustained load intensities
        
        W_q_vector_i = stats.gamma.ppf(stats.norm.cdf(W_q_vector_norm_i), k_q, scale=mew_U_q)

        # 4. Store this load case
       
        W_q_event_matrix[i, :] = W_q_vector_i
        
    # Generate intermittent load events
        
    N_p_events = np.random.poisson(v_p * T_yr)
   
    p_times = []
    t_current = 0

    while t_current < T_total:
        
          dt = np.random.exponential(scale=1 / v_p)
          t_current += dt
    
          if t_current < T_total:
            
             p_times.append(t_current)

    N_p_events = len(p_times)
    W_p_event_matrix = np.zeros((N_p_events, num_grid))

    for i in range(N_p_events):
        # 1. Generate independent standard normal vector
       
        W_p_uncorr_vector_i = np.random.normal(0, 1, size=num_grid)

        # 2. Apply transformation to get correlated standard normals
       
        W_p_vector_norm_i = L_p @ W_p_uncorr_vector_i

        # 3. Transform to gamma-distributed sustained load intensities
        
        W_p_vector_i = stats.gamma.ppf(stats.norm.cdf(W_p_vector_norm_i), k_p, scale=mew_U_p)

        # 4. Store this load case
       
        W_p_event_matrix[i, :] = W_p_vector_i
        
    # Assign each intermittent load vector to its correponding sustained load vector and sum
    
    p_by_q = {i: [] for i in range(len(starts_q))}

    for i, t_p in enumerate(p_times):
        
        for j, (start, end) in enumerate(zip(starts_q, ends_q)):
            
            if start <= t_p < end:
            
                p_by_q[j].append(i)
           
                break
    
    # Create a list for the total load vectors
    total_load_vectors = []

    # HERE IS WHERE YOU ADD THE NEW CODE:
    for j in range(len(durations_q)):
         W_q = W_q_event_matrix[j, :]  # sustained load vector for time segment j
    
         for i in p_by_q[j]:  # indices of intermittent events in this time window
            
             W_p = W_p_event_matrix[i, :]  # intermittent load vector
        
             # Calculate random factors for each grid point
             
             p_factors = []
        
             for grid_point in range(num_grid):
            
                 # Generate random location parameter
            
                 a_loc = np.random.normal(0.5, 0.5 * 0.2) * span
            
                 a_loc = np.clip(a_loc, 0, span)
            
                 b_loc = span - a_loc
            
                 # Calculate p_factor based on condition
            
                 if a_loc <= b_loc:
               
                    p_factor = 12 * a_loc * (b_loc ** 2) / (span ** 3)
            
                 else:  # a_loc > b_loc
                
                     p_factor = 12 * b_loc * (a_loc ** 2) / (span ** 3)
                
                 p_factors.append(p_factor)
        
             # Apply the random factors to W_p and combine with W_q
        
             W_p_adjusted = W_p * np.array(p_factors)
        
             W_total = W_q + W_p_adjusted
        
             # Store the total load vector
        
             total_load_vectors.append(list(W_total))
                 
    # Create a list within a list for the total load vectors
    
    List_of_total_load_cases = []
    
    for i in range(len(total_load_vectors)):
        
        List_of_total_load_on_each_grid = list(total_load_vectors[i])
        List_of_total_load_cases.append(List_of_total_load_on_each_grid)
            
    # Access a load to test
    
    for i in range(len(List_of_total_load_cases)):
        
        for j in range(num_grid):
            
            W_T = List_of_total_load_cases[i][j]
    
    
    
    # Loop through each grid point to define concrete random variables
    
    for j in range (num_grid): 
        
        # Extract correlated standard normal variable U for each point in the model
        
        U_ij = U_vector[j]
        
        # Generate basic concrete compressive strength f_co for each point
        
        # Take a single sample from Student's T distribution
        
        t_v = t.rvs(v)
        
        # Caclulate f_co
        
        f_co = np.exp(m + U_ij * sigma_j)
     
        # Calculate fc using JCSS PMC Eq.12, add a confinement factor to tune to abaqus model
        
        f_c_unconf = alpha * (f_co ** lambda_fc) * Y_1j
     
        conf_dev = 3.5 * ((conf_factor * f_c_unconf) ** (3/4)) * (f_c_unconf ** (1/4))
        
        f_c = f_c_unconf + conf_dev
       
        f_c_dict[j] = f_c
        
        # Calculate corresponding tensile strength according to JCSS PMC Eq.14
        
        f_ct = 0.3 * (f_c ** (2/3)) * Y_2j
        
        # Correlated steel properties
        
        # Generate reinforcement placement and properties for rib and scale and shift vectors
        
        A_s_list = []
        c_bot_list = []
        strain_s_y_list = []
        strain_s_ult_list = []
        bars = []
        
        # Cover for rib
        
        c_dev_bot = np.random.normal(c_bot_dev_mu, c_bot_dev_std)
        c_dev_top = np.random.normal(c_top_dev_mu, c_top_dev_std)
        
        c_bot = c + c_dev_bot
        c_top = c + c_dev_top
        
        c_bot_dict[j] = c_bot
        c_top_dict[j] = c_top
        
        if j % 2 != 0  and num_points <= j < num_points * 2: # mid span ribs are odd on second rung
        
           placement_opt = [num_bars_mid_p1, num_bars_mid_p2] # two different reinforcement possibilities for mid span
           
           placement_prob = [mid_p1, mid_p2] # probabilities for each scenario
        
           num_bars = np.random.choice(placement_opt, p = placement_prob) # choice based on probability
           
           if num_bars == num_bars_mid_p1:
               
               bar_diameter = bar_diameter_mid_p1
               
               for i in range(num_bars):
                   
                   A_s_nom = np.pi * ((bar_diameter / 2) ** 2) # calculate area of each bar
                   
                   # Young's Modulus of steel
                   
                   E_s = np.random.normal(E_s_mean, E_s_std)     # shows no correlation with other properties based on Wisniewski et al, 2009 Figure 7
                   
                   # generate vector of uncorrelated standard normals
                   
                   standard_steel_vector_uncorr = np.random.randn(4)
                   
                   # Apply matrix L to create correlated vector of steel properties
                   
                   standard_steel_vector = L_s @ standard_steel_vector_uncorr
                   
                   # Scale and shift to required distributions for each property
                   
                   A_s = standard_steel_vector[0] * A_s_std + A_s_nom
                   
                   f_y = standard_steel_vector[1] * f_y_std + f_y_mean
                  
                   f_ult = standard_steel_vector[2] * f_ult_std + f_ult_mean
                   
                   strain_s_ult = standard_steel_vector[3] * strain_s_ult_std + strain_s_ult_mean
                   
                   strain_s_y = f_y / E_s
                   
                   # Store properties in a dictionary
                   
                   bar = {
                       'A_s': A_s,
                       'f_y': f_y,
                       'f_ult': f_ult,
                       'E_s': E_s,
                       'strain_s_y': strain_s_y,
                       'strain_s_ult': strain_s_ult,
                       'bar_diameter': bar_diameter
                   }
                   bars.append(bar)
                   
                   A_s_list.append(A_s)
                   strain_s_y_list.append(strain_s_y)
                   strain_s_ult_list.append(strain_s_ult)
               
           else:
                
                bar_diameter = bar_diameter_mid_p2
                
                for i in range(num_bars):
                    
                       A_s_nom = np.pi * ((bar_diameter / 2) ** 2) # calculate area of each bar
                       
                       # Young's Modulus of steel
                       
                       E_s = np.random.normal(E_s_mean, E_s_std)     # shows no correlation with other properties based on Wisniewski et al, 2009 Figure 7
                       
                       # generate vector of uncorrelated standard normals
                       
                       standard_steel_vector_uncorr = np.random.randn(4)
                       
                       # Apply matrix L to create correlated vector of steel properties
                       
                       standard_steel_vector = L_s @ standard_steel_vector_uncorr
                       
                       # Scale and shift to required distributions for each property
                       
                       A_s = standard_steel_vector[0] * A_s_std + A_s_nom
                       
                       f_y = standard_steel_vector[1] * f_y_std + f_y_mean
                      
                       f_ult = standard_steel_vector[2] * f_ult_std + f_ult_mean
                       
                       strain_s_ult = standard_steel_vector[3] * strain_s_ult_std + strain_s_ult_mean
                       
                       strain_s_y = f_y / E_s
                  
                       bar = {
                           'A_s': A_s,
                           'f_y': f_y,
                           'f_ult': f_ult,
                           'E_s': E_s,
                           'strain_s_y': strain_s_y,
                           'strain_s_ult': strain_s_ult,
                           'bar_diameter': bar_diameter
                       }
                       bars.append(bar)
                       
                       A_s_list.append(A_s)
                       strain_s_y_list.append(strain_s_y)
                       strain_s_ult_list.append(strain_s_ult)
           
        if j % 2 == 0  and num_points <= j < num_points * 2: # slabs are even on second rung
        
            bar_diameter = bar_diameter_slab
            
            for i in range(num_bars_slab):
                
                A_s_nom = np.pi * ((bar_diameter / 2) ** 2) # calculate area of each bar
                
                # Young's Modulus of steel
                
                E_s = np.random.normal(E_s_mean, E_s_std)
                
                # generate vector of uncorrelated standard normals
                
                standard_steel_vector_uncorr = np.random.randn(4)
                
                # Apply matrix L to create correlated vector of steel properties
                
                standard_steel_vector = L_s @ standard_steel_vector_uncorr
                
                # Scale and shift to required distributions for each property
                
                A_s = standard_steel_vector[0] * A_s_std + A_s_nom
                
                f_y = standard_steel_vector[1] * f_y_std + f_y_mean
                
                f_ult = standard_steel_vector[2] * f_ult_std + f_ult_mean
                
                strain_s_ult = standard_steel_vector[3] * strain_s_ult_std + strain_s_ult_mean
                
                strain_s_y = f_y / E_s
                
                bar = {
                    'A_s': A_s,
                    'f_y': f_y,
                    'f_ult': f_ult,
                    'E_s': E_s,
                    'strain_s_y': strain_s_y,
                    'strain_s_ult': strain_s_ult,
                    'bar_diameter': bar_diameter
                }
                bars.append(bar)
                
                A_s_list.append(A_s)
                strain_s_y_list.append(strain_s_y)
                strain_s_ult_list.append(strain_s_ult)
                
        if j % 2 == 0 and j < num_points:# even ribs on first rung, consider end bar placement
            
            bar_diameter = [bar_diameter_end_p1, bar_diameter_end_p2]
            
            bar_diameter_prob = [end_p1, end_p2]
            
            bar_diameter = np.random.choice(bar_diameter, p = bar_diameter_prob)
            
            num_bars = num_bars_end
            
            for i in range(num_bars):
                
                A_s_nom = np.pi * ((bar_diameter / 2) ** 2) # calculate area of each bar
                
                # Young's Modulus of steel
                
                E_s = np.random.normal(E_s_mean, E_s_std)
                
                # generate vector of uncorrelated standard normals
                
                standard_steel_vector_uncorr = np.random.randn(4)
                
                # Apply matrix L to create correlated vector of steel properties
                
                standard_steel_vector = L_s @ standard_steel_vector_uncorr
                
                # Scale and shift to required distributions for each property
                
                A_s = standard_steel_vector[0] * A_s_std + A_s_nom
                
                f_y = standard_steel_vector[1] * f_y_std + f_y_mean
                
                f_ult = standard_steel_vector[2] * f_ult_std + f_ult_mean
                
                strain_s_ult = standard_steel_vector[3] * strain_s_ult_std + strain_s_ult_mean
                
                strain_s_y = f_y / E_s
                
                bar = {
                    'A_s': A_s,
                    'f_y': f_y,
                    'f_ult': f_ult,
                    'E_s': E_s,
                    'strain_s_y': strain_s_y,
                    'strain_s_ult': strain_s_ult,
                    'bar_diameter': bar_diameter
                }
                bars.append(bar)
                
                A_s_list.append(A_s)
                strain_s_y_list.append(strain_s_y)
                strain_s_ult_list.append(strain_s_ult)
            
        if j % 2 == 0 and num_points * 2 <= j: # last rung
            
            for i in range(num_bars):
                
                bar_diameter = [bar_diameter_end_p1, bar_diameter_end_p2]
                
                bar_diameter_prob = [end_p1, end_p2]
                
                bar_diameter = np.random.choice(bar_diameter, p = bar_diameter_prob)
                
                A_s_nom = np.pi * ((bar_diameter / 2) ** 2) # calculate area of each bar
                
                # Young's Modulus of steel
                
                E_s = np.random.normal(E_s_mean, E_s_std)
                
                # generate vector of uncorrelated standard normals
                
                standard_steel_vector_uncorr = np.random.randn(4)
                
                # Apply matrix L to create correlated vector of steel properties
                
                standard_steel_vector = L_s @ standard_steel_vector_uncorr
                
                # Scale and shift to required distributions for each property
                
                A_s = standard_steel_vector[0] * A_s_std + A_s_nom
                
                f_y = standard_steel_vector[1] * f_y_std + f_y_mean
                
                f_ult = standard_steel_vector[2] * f_ult_std + f_ult_mean
                
                strain_s_ult = standard_steel_vector[3] * strain_s_ult_std + strain_s_ult_mean
                
                strain_s_y = f_y / E_s
                
                bar = {
                    'A_s': A_s,
                    'f_y': f_y,
                    'f_ult': f_ult,
                    'E_s': E_s,
                    'strain_s_y': strain_s_y,
                    'strain_s_ult': strain_s_ult,
                    'bar_diameter': bar_diameter
                }
                bars.append(bar)
                
                A_s_list.append(A_s)
                strain_s_y_list.append(strain_s_y)
                strain_s_ult_list.append(strain_s_ult)
           
        if j % 2 != 0 and (j >= num_points * 2 or j < num_points): # end rungs slabs which aren't relevant
            
            A_s = 0
            f_y = 0
            f_ult = 0
            E_s = 0
            strain_s_y = 0
            strain_s_ult = 0
            bar_diameter = 0
            
            strain_s_y = 0
            
            bar = {
                'A_s': A_s,
                'f_y': f_y,
                'f_ult': f_ult,
                'E_s': E_s,
                'strain_s_y': strain_s_y,
                'strain_s_ult': strain_s_ult,
                'bar_diameter': bar_diameter
            }
            bars.append(bar)
            
            A_s_list.append(A_s)
            strain_s_y_list.append(strain_s_y)
            strain_s_ult_list.append(strain_s_ult)
            
        # Pick the relevant bar parameters for later calculations
                    
        strain_s_y_max = np.max(strain_s_y_list)
        strain_s_ult_dict.append(np.max(strain_s_ult_list))
        A_s_max = np.max(A_s_list)
        A_s_total = np.sum(A_s_list)
        A_s_total_dict[j] = np.sum(A_s_list)
        
        # Geometry - define slab depth for all grid points, but wont be needed for all makes it easier
                
        h_f_dev = np.random.normal(h_f_dev_mu, h_f_dev_std)
        
        h_f = h_f_nom + h_f_dev
        
        h_f_dict[j] = h_f
        
        # Mid span
        
        if j % 2 != 0  and num_points <= j < num_points * 2: # mid span ribs are odd on second rung
        
           # Dimensions
        
           h_dev = np.random.normal(h_dev_mu, h_dev_std)
        
           b_dev = np.random.normal(b_dev_mu, b_dev_std)
        
           h_f_dev = np.random.normal(h_f_dev_mu, h_f_dev_std)
        
           h = h_nom + h_dev
        
           b = b_nom + b_dev
        
           b_dict[j] = b
           h_dict[j] = h
           
        # End span 1
        
        if j % 2 == 0 and j < num_points:
            
            # Dimensions
         
            h_dev = np.random.normal(h_dev_mu, h_dev_std)
         
            b_dev = np.random.normal(b_dev_mu, b_dev_std)
         
            h_f_dev = np.random.normal(h_f_dev_mu, h_f_dev_std)
         
            h = h_nom + h_dev
         
            b = b_nom + b_dev
         
            b_dict[j] = b
            h_dict[j] = h
            
        # End span 2
        
        if j % 2 == 0 and j >= num_points * 2:
            
            # Dimensions
         
            h_dev = np.random.normal(h_dev_mu, h_dev_std)
         
            b_dev = np.random.normal(b_dev_mu, b_dev_std)
         
            h_f_dev = np.random.normal(h_f_dev_mu, h_f_dev_std)
         
            h = h_nom + h_dev
         
            b = b_nom + b_dev
         
            b_dict[j] = b
            h_dict[j] = h
            
        # Effective section size for mid ribs
        
        if j % 2 != 0  and num_points <= j < num_points * 2: # mid span ribs are odd on second rung
        
           b = b_dict[j]
           
           l_0 = l_0_factor * span
           
           b_eff = b + 2 * min(0.2 * ((rib_spacing - b) / 2) + 0.1 * l_0, 0.2 * l_0, (rib_spacing - b)/ 2)  # equation for effective T section from Eurocode 2
        
           b_eff_dict[j] = b_eff
           
        # End span 1
        
        if j % 2 == 0 and j < num_points:
            
            b = b_dict[j]
            
            b_eff = rib_spacing
            
            b_eff_dict[j] = b_eff
                      
        # End span 2
        
        if j % 2 == 0 and j >= num_points * 2:
            
            b = b_dict[j]
            
            b_eff = rib_spacing
            
            b_eff_dict[j] = b_eff
            
    
        # Concrete Elastic Modulus
        
        beta_d = 0.7
    
        # Creep coefficient Eurocode 2 appendix B
    
        # Basic creep coefficient
      
        beta_bc_fc = 1.8 / (f_c ** 0.7) # accounts for effect of concrete strength on basic creep
 
        t_load_adj = t_load * (9 / (2 + (t_load ** 1.2) + 1) ** alpha_sc)
    
        if t_load_adj <= 0.5:
        
           t_load_adj = 0.5
    
        beta_bc_t = np.log((((30 / t_load_adj) + 0.35) ** 2) * (tau - t_load) + 1) # time development of basic creep
        
        creep_basic = beta_bc_fc * beta_bc_t
      
        # Drying creep coefficient
    
        beta_dc_fc = 412 / (f_c ** 1.4) # accounts for effect of concrete strength on drying creep
        
        if j % 2 == 0:
        
           A_c = h_f * b_eff + h * b # concrete cross-sectional area
        
           u = b_eff + 2 * h + (b_eff - b) + b # perimeter in contact with atmosphere (everything except the sides of the flange)
           
        else:
            
            A_c = span * h_f # concrete cross-sectional area
            
            u = 2 * (rib_spacing - b)  # perimeter in contact with atmoshpere (everything except the sides of the flange)
            
        h_n = 2 * A_c / u # notional size mm
            
        beta_dc_RH = (1 - (RH / 100)) / ((0.1 * (h_n / 100)) ** (1/3)) # with RH already in decimal form accounts for the effect of RH on drying creep
        
        beta_dc_t = 1 / (0.1 + (t_load_adj ** 0.2)) 
        
        gamma_t = 1 / (2.3 + (3.5 / (t_load_adj ** (1/2))))
        
        alpha_fc = (35 / f_c) ** 0.5
        
        beta_h = 1.5 * h_n + 250 * alpha_fc
        
        if beta_h <= 1500 * alpha_fc:
            
            beta_h = 1500 * alpha_fc
            
        beta_dc_t_tau = ((tau - t_load) / (beta_h + (tau - t_load))) ** gamma_t
        
        creep_drying = beta_dc_fc * beta_dc_RH * beta_dc_t * beta_dc_t_tau
      
        # Calculate creep factor from basic and drying components
        
        creep_factor = creep_basic + creep_drying
       
        creep_factor_dict[j] = creep_factor
        
        # Work out concrete Young's Modulus according to JCSS PMC Section 3
        
        E_c = 10.5 * (f_c ** (1/3)) * Y_3j * (1 / (1 + beta_d * creep_factor)) * 1000 # convert to MPa
      
        E_c_dict[j] = E_c
        
        # Concrete ultimate strain calculation
        
        strain_c_ult = 6 * (10 ** (-3)) * (f_c ** (-1/6)) * Y_4j * (1 + beta_d * creep_factor)
        
        strain_c_ult_dict[j] = strain_c_ult
        
        # Determine if rib is brittle or ductile 
        
        if j % 2 != 0  and num_points <= j < num_points * 2: # consider mid span ribs only
        
           # Assuming rectangular compressive stress block for simplicity Eurocode 2 8.1.2

           # For a balanced section strains are balanced ultimate concrete strain = yield steel strain
        
           # Consider bar (if multiple bars) with greatest yield strain
           
           # Strain compatibility
           
           bar_diameter_max = 2 * np.sqrt(A_s_max / np.pi)

           d = h + h_f - c_bot - (bar_diameter_max / 2) # effective depth from top of the T section
           
           d_dict[j] = d
           
           # Transformed section method for cracked elastic section
           
           n = E_s / E_c # modular ratio calculation
           
           A_s_eff = n * A_s_total # effective area of steel
           
           # Reference at top of section for first moment of area centroid

           a_quad = b_eff / 2
           
           b_quad = A_s_eff
           
           c_quad = - A_s_eff * d
           
           disc = (b_quad ** 2) - 4 * a_quad * c_quad
           
           x_sls = (-b_quad + np.sqrt(disc)) / (2 * a_quad) 
         
           # Check for NA in the web 
           
           if x_sls > h_f:
               
              # Need to recheck for some area in web
              
              A_c_f = b_eff * h_f
              
              a_quad = b
              
              b_quad = 2 * (h_f * (b_eff - b) + A_s_eff)
              
              c_quad  = - ((h_f ** 2) * (b_eff - b) + 2 * A_s_eff * d)
              
              disc = (b_quad ** 2) - 4 * a_quad * c_quad
              
              disc = (b_quad ** 2) - 4 * a_quad * c_quad
              
              x_sls = (-b_quad + np.sqrt(disc)) / (2 * a_quad) 
             
            # Find transformed second moment of area about x
            
           if x_sls <= h_f:
               
               y_conc = x_sls / 2 # half depth of compressive section of flange
               
               y_steel = d  - x_sls
               
               I_eff = (b_eff * (x_sls ** 3) / 12) + b_eff * x_sls * ((x_sls - y_conc) ** 2) + A_s_eff * ((x_sls - y_steel) ** 2)
           
           if x_sls > h_f:
               
               y_f = h_f / 2
               
               y_w = h_f + (h_f - x_sls) / 2
               
               y_steel = d - x_sls
               
               I_eff = (b_eff * (h_f ** 3) / 12) + b_eff * h_f * ((x_sls - y_f) ** 2) + (b * (((x_sls - h_f) / 2) ** 3) / 12) + b * (x_sls - h_f) / 2 * ((x_sls - y_w) ** 2) + A_s_eff * ((x_sls - y_steel) ** 2)
               
           # Calculate effective stiffness of transformed section at mid-span- only accurate for serviceability state where load is low
           
           EI = E_c * I_eff
           
           EI_dict[j] = EI

           
           # Use strain compatibility to find neutral axis depth of balanced section
           
           x_balanced = (d * strain_c_ult) / (strain_c_ult + strain_s_y_max) # According to ACI 318-99
           
           if x_balanced <= x_sls: # section is balanced or over-reinforced so brittle
           
              sec_brittle_dict[j] = True
              
           # Now consider ultimate limit state assuming rectangular stress block - calculate the ultimate moment assuming whitney stress block for ductile section
           
           if sec_brittle_dict[j] == False: # tension controlled section
           
              # Assume steel has fully yielded and concrete is at point of yielding
           
              # First assume neutral axis is in the flange at this point
           
              x_ult = (sum(bar['f_y'] * bar['A_s'] for bar in bars)) / (0.8 * 0.85 * f_c * b_eff)
           
              if x_ult > h_f:
                  
                  # Repeat calculation considering part of the web in compression
                  
                  F_c_flange = b_eff * h_f * 0.85 * f_c # compressive force carried in the flange
                  
                  F_t = sum(bar['f_y'] * bar['A_s'] for bar in bars) # tensile force in the rebar
                  
                  F_c_web = F_t - F_c_flange # compressive force required in the web
                  
                  x_web = F_c_web / (0.85 * f_c * b)
                 
                  x_ult = (h_f + x_web) * 1.25 # 0.8x stress block adjustment
                  
              z_ult = d - 0.4 * x_ult
              
              z_dict[j] = z_ult
                
              M_ult_capacity = z_ult * (sum(bar['f_y'] * bar['A_s'] for bar in bars))
              
              M_capacity_ult_dict[j] = M_ult_capacity
              
              M_capacity = (d - x_sls * 0.5) * (sum(bar['f_y'] * bar['A_s'] for bar in bars))
              
              M_capacity_dict[j] = M_capacity
              
           if sec_brittle_dict[j] == True: # compression controlled beam
           
               # Find ultimate capacity of the compression controlled beam
               
               # Assume whiteny stress block - but neutral axis at serviceability- store this as just capacity
               
               M_capacity = 0.6 * x_sls * 0.85 * f_c
               
               M_capacity_dict[j] = M_capacity
               
           # Calculate vertical shear capacity, note that shear failure should be treated as brittle in the model even if the rib is classified as ductile for flexure
           
           # Minimum shear stress resistance calculation according to Eurocode 2
           
           # Omitting partial safety factors as this is a probabilistic approach
           
           # Note T_section has geometric discontinutities
           
           # Use simple conservative approach for calculating the minimum vertical shear capacity
           
           D_lower = 16 # smallest value of upper sieve size D, assume 20 for now
           
           d_dg = 16 + D_lower # size parameter for failure zone roughness
           
           if d_dg >= 40:
               
               d_dg = 40
               
           v_capacity = 11 * np.sqrt((f_c / f_y) * (d_dg / d))
           
           v_capacity_dict[j] = v_capacity
           
           # Load capacity for flexure
           
           node_load_max_bending = (24 * M_capacity) / (span) # assuming support condition for encastre
       
           # Load capacity for vertical shear
           
           node_load_max_shear = v_capacity * 0.9 * b * z_ult
           
           node_load_max = min(node_load_max_bending, node_load_max_shear)
           
           node_load_max_dict[j] = node_load_max
        
        if j % 2 == 0  and num_points <= j < num_points * 2: # slab
         
            # Calculate slab second moment of area assuming it remains in serviceability limit state
            
            # Transformed section method for cracked elastic section
            
            d = h_f - c_bot / 2
            
            n = E_s / E_c # modular ratio calculation
            
            A_s_eff = n * A_s_total # effective area of steel
            
            # Reference at top of section for first moment of area centroid

            a_quad = rib_spacing / 2
            
            b_quad = A_s_eff
            
            c_quad = - A_s_eff * d
            
            disc = (b_quad ** 2) - 4 * a_quad * c_quad
            
            x_sls = (-b_quad + np.sqrt(disc)) / (2 * a_quad) 
            
            y_conc = x_sls / 2 # half depth of compressive section of flange
            
            y_steel = d  - x_sls
            
            I_eff = (span * (x_sls ** 3) / 12) + span * x_sls * ((x_sls - y_conc) ** 2) + A_s_eff * ((x_sls - y_steel) ** 2)
         
            I_dict[j] = I_eff
            
        if j % 2 == 0 and j < num_points: # end span 1
        
           # Transformed section method for cracked elastic section
           
           d = h + h_f - c_top / 2
           
           n = E_s / E_c # modular ratio calculation
           
           A_s_eff = n * A_s_total # effective area of steel
           
           # Reference at top of section for first moment of area centroid

           a_quad = rib_spacing / 2
           
           b_quad = A_s_eff
           
           c_quad = - A_s_eff * d
           
           disc = (b_quad ** 2) - 4 * a_quad * c_quad
           
           x_sls = (-b_quad + np.sqrt(disc)) / (2 * a_quad) 
           
           y_conc = x_sls / 2 # half depth of compressive section of flange
           
           y_steel = d  - x_sls
           
           I_eff = (b_eff * (x_sls ** 3) / 12) + b_eff * x_sls * ((x_sls - y_conc) ** 2) + A_s_eff * ((x_sls - y_steel) ** 2) 
           
           # Calculate effective stiffness of transformed section at mid-span- only accurate for serviceability state where load is low
           
           EI = E_c * I_eff 
           
           EI_dict[j] = EI
           
           # Use strain compatibility to find neutral axis depth of balanced section
           
           x_balanced = (d * strain_c_ult) / (strain_c_ult + strain_s_y_max) # According to ACI 318-99
           
           if x_balanced <= x_sls: # section is balanced or over-reinforced so brittle
           
              sec_brittle_dict[j] = True
              
           # Now consider ultimate limit state assuming rectangular stress block - calculate the ultimate moment assuming whitney stress block for ductile section
           
           if sec_brittle_dict[j] == False: # tension controlled section
           
              # Assume steel has fully yielded and concrete is at point of yielding
           
              # First assume neutral axis is in the flange at this point
           
              x_ult = (sum(bar['f_y'] * bar['A_s'] for bar in bars)) / (0.8 * 0.85 * f_c * b_eff)
                  
              z_ult = d - 0.4 * x_ult
              
              z_dict[j] = z_ult
                
              M_ult_capacity = z_ult * (sum(bar['f_y'] * bar['A_s'] for bar in bars))
              
              M_capacity_ult_dict[j] = M_ult_capacity
              
              M_capacity = (d - x_sls * 0.5) * (sum(bar['f_y'] * bar['A_s'] for bar in bars))
              
              M_capacity_dict[j] = M_capacity
           
           if sec_brittle_dict[j] == True: # compression controlled beam
           
               # Find ultimate capacity of the compression controlled beam
               
               # Assume whiteny stress block - but neutral axis at serviceability
               
               M_capacity = 0.6 * x_sls * 0.85 * f_c
               
               M_capacity_dict[j] == M_capacity
               
           # Calculate vertical shear capacity, note that shear failure should be treated as brittle in the model even if the rib is classified as ductile for flexure
           
           # Minimum shear stress resistance calculation according to Eurocode 2
           
           # Omitting partial safety factors as this is a probabilistic approach
           
           # Note T_section has geometric discontinutities
           
           # Use simple conservative approach for calculating the minimum vertical shear capacity
           
           D_lower = 20 # smallest value of upper sieve size D, assume 20 for now
           
           d_dg = 16 + D_lower # size parameter for failure zone roughness
           
           if d_dg >= 40:
               
               d_dg = 40
               
           v_capacity = 11 * np.sqrt((f_c / f_y) * (d_dg / d))
           
           v_capacity_dict[j] = v_capacity
           
           # Load capacity for flexure
           
           node_load_max_bending = (24 * M_capacity) / (span) # assuming support condition for encastre
         
           # Load capacity for vertical shear
           
           node_load_max_shear = 2 * v_capacity * 0.9 * b * z_ult
           
           node_load_max = min(node_load_max_bending, node_load_max_shear)
           
           node_load_max_dict[j] = node_load_max
       
        if j % 2 == 0 and j >= num_points * 2: # end span 2
        
            # Transformed section method for cracked elastic section
            
            d = h + h_f - c_top / 2
            
            n = E_s / E_c # modular ratio calculation
            
            A_s_eff = n * A_s_total # effective area of steel
            
            # Reference at top of section for first moment of area centroid

            a_quad = rib_spacing / 2
            
            b_quad = A_s_eff
            
            c_quad = - A_s_eff * d
            
            disc = (b_quad ** 2) - 4 * a_quad * c_quad
            
            x_sls = (-b_quad + np.sqrt(disc)) / (2 * a_quad) 
            
            y_conc = x_sls / 2 # half depth of compressive section of flange
            
            y_steel = d  - x_sls
            
            I_eff = (b_eff * (x_sls ** 3) / 12) + b_eff * x_sls * ((x_sls - y_conc) ** 2) + A_s_eff * ((x_sls - y_steel) ** 2) 
             
            # Calculate effective stiffness of transformed section at mid-span- only accurate for serviceability state where load is low
            
            EI = E_c * I_eff 
            
            EI_dict[j] = EI
          
            # Use strain compatibility to find neutral axis depth of balanced section
            
            x_balanced = (d * strain_c_ult) / (strain_c_ult + strain_s_y_max) # According to ACI 318-99
            
            if x_balanced <= x_sls: # section is balanced or over-reinforced so brittle
            
               sec_brittle_dict[j] = True
               
            # Now consider ultimate limit state assuming rectangular stress block - calculate the ultimate moment assuming whitney stress block for ductile section
            
            if sec_brittle_dict[j] == False: # tension controlled section
            
               # Assume steel has fully yielded and concrete is at point of yielding
            
               # First assume neutral axis is in the flange at this point
            
               x_ult = (sum(bar['f_y'] * bar['A_s'] for bar in bars)) / (0.8 * 0.85 * f_c * b_eff)
                   
               z_ult = d - 0.4 * x_ult
               
               z_dict[j] = z_ult
                 
               M_ult_capacity = z_ult * (sum(bar['f_y'] * bar['A_s'] for bar in bars))
               
               M_capacity_ult_dict[j] = M_ult_capacity
               
               M_capacity = (d - x_sls * 0.5) * (sum(bar['f_y'] * bar['A_s'] for bar in bars))
               
               M_capacity_dict[j] = M_capacity
               
            if sec_brittle_dict[i] == True: # compression controlled beam
            
                # Find ultimate capacity of the compression controlled beam
                
                # Assume whiteny stress block - but neutral axis at serviceability
                
                M_capacity = 0.6 * x_sls * 0.85 * f_c
                
                M_capacity_dict[j] = M_capacity
          
            # Calculate vertical shear capacity, note that shear failure should be treated as brittle in the model even if the rib is classified as ductile for flexure
            
            # Minimum shear stress resistance calculation according to Eurocode 2
            
            # Omitting partial safety factors as this is a probabilistic approach
            
            # Note T_section has geometric discontinutities
            
            # Use simple conservative approach for calculating the minimum vertical shear capacity
            
            D_lower = 20 # smallest value of upper sieve size D, assume 20 for now
            
            d_dg = 16 + D_lower # size parameter for failure zone roughness
            
            if d_dg >= 40:
                
                d_dg = 40
                
            v_capacity = 11 * np.sqrt((f_c / f_y) * (d_dg / d))
            
            v_capacity_dict[j] = v_capacity
            
            # Load capacity for flexure
            
            node_load_max_bending = (12 * M_capacity) / (span) # assuming support condition for encastre
          
            # Load capacity for vertical shear
            
            node_load_max_shear = 2 * (v_capacity * 0.9 * b * z_ult) 
            
            node_load_max = min(node_load_max_bending, node_load_max_shear)
           
            node_load_max_dict[j] = node_load_max
            
    # Calculate equivalent translational stiffness
    
    for j in range(num_grid):
        
        if j % 2 != 0  and num_points <= j < num_points * 2: # consider mid span ribs first
        
           EI_1 = EI_dict[j - num_points]
           EI_2 = EI_dict[j]
           EI_3 = EI_dict[j + num_points]
           
           EI_eq = span / ((pot / EI_1) + ((span - 2 * pot) / EI_2) + (pot / EI_3))
           
           # Calculate equiavlent translational stiffness for mid span assuming encastre boundary conditions
           
           k = 384 * EI_eq / (span ** 3) # equation for mid span deflection with the total UDL applied, which is equivalent to the real scenario rearranged
           
           k_dict[j] = k 
           
           EI_eq_dict[j] = EI_eq
           
           original_EI[j] = EI_eq
           
      # Calculate total load for each of the load cases
             
      # Calculate self-weight on each rib as a nodel load in N
      
      # Calculate slab self weight load on each rib
         
    for j in range(num_grid):
           
           f_c = f_c_dict[j]
           
           b = b_dict[j]
           
           h = h_dict[j]
           
           h_f = h_f_dict[j]
           
           b_eff = b_eff_dict[j]
           
           mass_den_conc = 8.602 * f_c + 2110 # relationship derived by Iffat 2015 for curshed stone aggregate concrete in an appropriate strength range for this application
           
           den_conc = mass_den_conc * 9.81e-9 # N/mm^3
      
           if j % 2 != 0  and num_points <= j < num_points * 2: # consider mid span ribs only
          
              sw_rib_node = den_conc * (b * h + b_eff * h_f) * span
              
              sw_node_load_dict[j] = sw_rib_node
              
           # Calculate nodal dead load
            
           d_load_node = d_load * rib_spacing * span
            
           # Create a new list to store the total loads with self-weight and dead load
        
           Total_node_load_cases = []

           # Add the selfweight and dead load to each total load value
        
           for i in range(len(List_of_total_load_cases)):
        
              # Create a new list for this load case
      
              node_loads_for_this_case = []
      
              for j in range(num_grid):
          
               # Calculate total node load for this grid point
          
                   node_load = (List_of_total_load_cases[i][j] / 1000) * rib_spacing * span + d_load_node + sw_node_load_dict[j]
          
                   # Add to the list for this case
          
                   node_loads_for_this_case.append(node_load)
              
              # Add the completed list for this case to the main list
         
              Total_node_load_cases.append(node_loads_for_this_case)
            
    for j in range(num_grid):
        
            # Compare lists of load and max load (before redistribution is considered) to find which rib will likely fail first
           
            if j % 2 != 0  and num_points <= j < num_points * 2: # consider mid span ribs first
            
               # Longitudinal shear capacity 
            
               # Transverse shear reinforcement
               
               if j == num_points:
                
                  A_s_total = A_s_total_dict[j+1] # just consider slab to the right
                
               if j == 2 * num_points - 1: 
                
                  A_s_total = A_s_total_dict[j-1] # just consider slab to the left
                
               if num_points < j < 2 * num_points - 1:
                
                  A_s_total = (A_s_total_dict[j+1] +  A_s_total_dict[j-1]) / 2 # average area across the flange
                
               # Calculate steel area - spacing ratio
              
               tau_capacity = A_s_total / spacing_slab
            
               tau_capacity_dict[j] = tau_capacity
               
    min_scale_factor = float('inf')  # Start with infinity



    # Check if Total_node_load_cases is empty
    if len(Total_node_load_cases) == 0:
        print("Total_node_load_cases is empty. Skipping to next simulation.")
        # Skip directly to the next simulation iteration
        # You might need to adjust these values depending on your simulation structure
        System_ULS_fail_counter += 0  # Don't count as failure
        System_SLS_fail_counter += 0  # Don't count as failure
    else:
        # There are load cases to process, proceed normally
        critical_load_case = 0 

        # Loop through each load case
        for i in range(len(Total_node_load_cases)):
               
            # Track the minimum scale factor for this specific load case
            min_scale_for_this_case = 1
            
            # Loop through each grid point
            for j in range(num_grid):
                        
                if j % 2 != 0 and num_points <= j < num_points * 2:
                           
                    # Calculate actual load at this point for this case
                    actual_load = Total_node_load_cases[i][j]
                    
                    # Skip division by zero
                    if actual_load == 0:
                        continue
                    
                    # Calculate the scale factor (ratio of max allowable to actual)
                    # Note: A scale factor < 1 means the actual load exceeds the maximum allowable
                    scale_factor_mid = node_load_max_dict[j] / actual_load 
                    scale_factor_end_1 = node_load_max_dict[j - num_points] / actual_load 
                    scale_factor_end_2 = node_load_max_dict[j + num_points] / actual_load 
                   
                    scale_factor = min(scale_factor_mid, scale_factor_end_1, scale_factor_end_2)
           
                    # Update the minimum scale factor for this load case if needed
                    if scale_factor < min_scale_for_this_case:
                        print(f"New critical case found: {i} with scale factor {scale_factor}")
                        min_scale_for_this_case = scale_factor

            # If this load case has a lower minimum scale factor than what we've seen before,
            # it becomes our critical case
            if min_scale_for_this_case < min_scale_factor:
                min_scale_factor = min_scale_for_this_case
                critical_load_case = i
                       
        print('CLC', critical_load_case)
                       
        scale_factor = min_scale_factor
                
        if scale_factor >= 1:
            scale_factor = 1
        
        # Adjust all nodal loads by scale factor
        node_load_dict = list(Total_node_load_cases[critical_load_case])
                 
        scaled_node_load_dict = [scale_factor * node_load for node_load in node_load_dict]

        print("Critical load case:", critical_load_case)
            
        # First branch: if scale factor is less than one 
        if scale_factor < 0.999:
            System_ULS_fail_counter, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter = iterate_full_load(System_ULS_fail_counter, k_dict, rib_ductile_fail_dict, rib_brittle_fail_dict, rib_status_dict, E_c_dict, I_dict, scaled_node_load_dict, strain_c_ult_dict, Lshear_counter, Vshear_counter, Flexure_fail_counter_mid, Flexure_fail_counter_end, collapse_counter)
                     
        # Check serviceability
        deflection_max = max(deflection_dict)
          
        print('max deflection = ', deflection_max)
          
        deflection_limit = span / 250
          
        print('deflection_limit = ', deflection_limit)
            
        if deflection_max > deflection_limit:
            System_SLS_fail_counter += 1
            
            
 # Continue with the rest of your code outside the if/else block
 # System reliability
print(System_ULS_fail_counter)

 # Calculate system serviceability reliability
serviceability_pf = System_SLS_fail_counter / num_sims
serviceability_reliability = -norm.ppf(serviceability_pf)
print('Probability of system serviceability failure = ', serviceability_pf)
print('Serviceability system reliability = ', serviceability_reliability)

 # Calculate system ultimate reliability
ultimate_pf = System_ULS_fail_counter / num_sims
ultimate_reliability = -norm.ppf(ultimate_pf)
print('Probability of system ultimate failure = ', ultimate_pf)
print('Ultimate system reliability = ', ultimate_reliability) 
print('Bending failures at end span = ', Flexure_fail_counter_end)
print('Bending failures at mid span = ', Flexure_fail_counter_mid) 
print('V. shear failures at end span = ', Vshear_counter)   
print('L. shear failures at mid span = ', Lshear_counter)  
print('Collapse of system failures = ', collapse_counter)     
   
 # End timing
end_time = time.time()
 # Calculate and print execution time
execution_time = end_time - start_time
print(f"Execution time: {execution_time:.4f} seconds")

print("Failed Rib Info (rib index, original EI):")
for rib_index, ei_value in failed_rib_info:
    print(f"Rib {rib_index}: Original EI = {ei_value:.2f}")
          
# Check each load against each rib capacity
    
# Initialize variables to track the minimum scale factor and corresponding load case
         
