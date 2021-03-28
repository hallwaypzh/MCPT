#include "TriangleMeshes.h"

#include <iostream>
#include <sstream>

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

TriangleMeshes::TriangleMeshes(std::string &dir, std::string &file_name) {
    std::string obj_path = dir + file_name;
    std::string mtl_path = dir;
    tinyobj::ObjReaderConfig reader_config;
    reader_config.mtl_search_path = dir;
    tinyobj::ObjReader reader;
    if (!reader.ParseFromFile(obj_path, reader_config)) {
        if (!reader.Error().empty()) {
            std::cerr << "TinyObjReader: " << reader.Error();
        }
        exit(1);
    }
    if (!reader.Warning().empty()) {
        std::cout << "TinyObjReader: " << reader.Warning();
    }
    auto& attrib = reader.GetAttrib();
    auto& shapes = reader.GetShapes();
    auto& materials = reader.GetMaterials();
    auto& vertices = attrib.vertices;
    auto& normals = attrib.normals;
    auto& texcoord = attrib.texcoords;  

    for (auto &m : materials) {
        std::string name = m.name;
        std::string diffuse_texture_name = m.diffuse_texname;
        glm::vec3 ka = glm::vec3(m.ambient[0], m.ambient[1], m.ambient[2]);
        glm::vec3 kd = glm::vec3(m.diffuse[0], m.diffuse[1], m.diffuse[2]);
        glm::vec3 ks = glm::vec3(m.specular[0], m.specular[1], m.specular[2]);
        Material material(name, ka, kd, ks, m.shininess, m.ior, m.illum);
        
        //load diffuse texture
        if (m.diffuse_texname.size() > 0) {
            //load diffuse texture map
            Texture tex(dir, diffuse_texture_name);
            material.diffuse_tex = tex;
        }
        
        if (m.unknown_parameter.size() > 0) {
            auto it = m.unknown_parameter.find("Le");
            if (it != m.unknown_parameter.end()) {
                //light source
                glm::vec3 le;
                std::stringstream ss((*it).second );
                for (int k = 0; k < 3; k++) {
                    ss >> le[k];
                    le[k] *= 1.f;
                }
                material.setLe(le);
            } else {
                std::cerr << "unkown parameters: " << std::endl;
            }
        }
        this->materials.push_back(material);
    }

    for (int i = 0; 3 * i < vertices.size(); i++) {
        this->vertices.push_back(glm::vec4(vertices[3 * i + 0], vertices[3 * i + 1], vertices[3 * i + 2], 1.f));
    }
    
    for (int i = 0; 3 * i < normals.size(); i++) {
        this->normals.push_back(glm::vec4(normals[3 * i + 0], normals[3 * i + 1], normals[3 * i + 2], 0.f));
    }

    for (int i = 0; 2 * i < texcoord.size(); i++) {
        this->texcoord.push_back(glm::vec2(texcoord[2 * i + 0], texcoord[2 * i + 1]));
    }

    // Loop over shapes
    
    for (size_t s = 0; s < shapes.size(); s++) {
        // Loop over faces(polygon)

        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            int fv = shapes[s].mesh.num_face_vertices[f];
            tinyobj::index_t id0 = shapes[s].mesh.indices[index_offset];
            // Loop over vertices in the face.
            for (size_t v = 1; v < fv - 1; v++) {
                // access to vertex
                tinyobj::index_t idi = shapes[s].mesh.indices[index_offset + v];
                tinyobj::index_t idj = shapes[s].mesh.indices[index_offset + v + 1];
                vids.push_back(std::make_tuple(id0.vertex_index, idi.vertex_index, idj.vertex_index));
                tids.push_back(std::make_tuple(id0.texcoord_index, idi.texcoord_index, idj.texcoord_index));
                nids.push_back(std::make_tuple(id0.normal_index, idi.normal_index, idj.normal_index));
            }
            mids.push_back(shapes[s].mesh.material_ids[f]);
            index_offset += fv;
        }
    }
    std::cout << "loaded scene successfully:)\n";
    std::cout << "num vertices: " << this->vertices.size() << std::endl;;
    std::cout << "num triangles: " << this->vids.size() << std::endl;
    std::cout << "num materials: " << this->materials.size() << std::endl;
}

void TriangleMeshes::applyTransform(glm::mat4x4 &modelview) {
    //transform vertices 
    for (auto &v : vertices) {
        v = glm::vec3(modelview * glm::vec4(v, 1.f));
    }

    // TO-DO transform normals
    // since modelview  is just view matrix, which is orthographic,
    // so we actually don't have to do anything
    glm::mat3x3 M(modelview);
    glm::mat3x3 invMT = glm::transpose(glm::inverse(M));
    for (auto &n : normals) {
        n = invMT * n;
    }
}

void TriangleMeshes::initializeTriangles(std::vector<Triangle> &triangles) {
    for (int i = 0; i < vids.size(); i++) {
        std::tuple<int, int, int> vid = vids[i];
        std::tuple<int, int, int> tid = tids[i];
        std::tuple<int, int, int> nid = nids[i];
        /*
        std::cout << m_vertices[std::get<0>(vids)].x << " " << m_vertices[std::get<0>(vids)].y << " " << m_vertices[std::get<0>(vids)].z << std::endl;
        std::cout << m_vertices[std::get<1>(vids)].x << " " << m_vertices[std::get<1>(vids)].y << " " << m_vertices[std::get<1>(vids)].z << std::endl;
        std::cout << m_vertices[std::get<2>(vids)].x << " " << m_vertices[std::get<2>(vids)].y << " " << m_vertices[std::get<2>(vids)].z << std::endl;
        */
        
        Triangle t(vertices[std::get<0>(vid)], vertices[std::get<1>(vid)], vertices[std::get<2>(vid)],
                   normals[std::get<0>(nid)], normals[std::get<1>(nid)], normals[std::get<2>(nid)],
                   texcoord[std::get<0>(tid)], texcoord[std::get<1>(tid)], texcoord[std::get<2>(tid)],
                   &materials[mids[i]]);
        if (materials[mids[i]].is_light)
            this->luminaires.push_back(t);
        triangles.push_back(t);
    }
}